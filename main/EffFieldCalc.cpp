/*
    tetmag - A general-purpose finite-element micromagnetic simulation software package
    Copyright (C) 2016-2023 CNRS and Université de Strasbourg

    Author: Riccardo Hertel

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU Affero General Public License as
    published by the Free Software Foundation, either version 3 of the
    License, or (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Affero General Public License for more details. 
 
    Contact: Riccardo Hertel, IPCMS Strasbourg, 23 rue du Loess, 
    	     67034 Strasbourg, France.
	     riccardo.hertel@ipcms.unistra.fr
	     
    You should have received a copy of the GNU Affero General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/
    
 /*
 * EffFieldCalc.cpp
 *
 *  Created on: May 25, 2021
 *      Author: riccardo
 */

#include <EffFieldCalc.h>
#include <iostream>
#include "auxiliaries.h"
using namespace Eigen;

enum Coords { x,	y,	z };

EffFieldCalc::EffFieldCalc(SimulationData &sd, const MeshData &msh) : nx(sd.nx), XC_field_OP(msh.stiff),
																	  A(sd.A), D(sd.D), Ku1(sd.Ku1), KuAxis(sd.Kuni),
																	  cubicAxes(sd.cubicAxes), Kc1(sd.Kc1), Kc2(sd.Kc2), Ksn(sd.Ks),
																	  nodeSurfaceArea(msh.nodeArea), nv_nx(msh.nv_nx),
																	  tGradX(-msh.tGradX), tGradY(-msh.tGradY), tGradZ(-msh.tGradZ),
																	  gradX(msh.gradX), gradY(msh.gradY),gradZ(msh.gradZ)
{
	Hexc = MatrixXd::Zero(nx, 3);
	Hani = MatrixXd::Zero(nx, 3);
	Hsrf = MatrixXd::Zero(nx, 3);
	Hdmi = MatrixXd::Zero(nx, 3);
	Hcub = MatrixXd::Zero(nx, 3);

	NodeVol = msh.NodeVolume;
	invNodeVol = NodeVol.cwiseInverse();
	setExchangeFieldOperator();
	setUniaxialAnisotropy();
	setSurfaceAnisotropy();
}

void EffFieldCalc::setExchangeFieldOperator()
{
	XC_field_OP = (A.cwiseProduct(invNodeVol)).asDiagonal() * XC_field_OP;
	XC_field_OP = (-2.) * XC_field_OP;
}

double EffFieldCalc::calcExchangeEnergy(const MRef &Mag)
{
	Hexc = exchangeField(Mag);
	return (NodeVol.transpose() * Hexc.cwiseProduct(Mag)).sum() / (-2.);
}

MatrixXd EffFieldCalc::uniaxialAnisotropyField(const MRef &Mag)
{
#ifdef _OPENMP
#pragma omp parallel for
	for (int i = 0; i < nx; ++i)
	{
		const Vector3d axis = KuAxis.row(i);
		Hani.row(i) = 2. * Ku1(i) * axis.dot(Mag.row(i)) * axis;
	}
#else
	Hani = 2. * Ku1.cwiseProduct(KuAxis.cwiseProduct(Mag).rowwise().sum()).asDiagonal() * KuAxis;
#endif
	return Hani;
}

void EffFieldCalc::setUniaxialAnisotropy()
{
	uniaxialAnisotropyEnergyOffset = NodeVol.cwiseProduct(Ku1).sum();
}

double EffFieldCalc::uniaxialEnergy(const MRef &Mag)
{
	VectorXd t1 = KuAxis.cwiseProduct(Mag).rowwise().sum();
	return uniaxialAnisotropyEnergyOffset - (Ku1.cwiseProduct(NodeVol)).cwiseProduct(t1).dot(t1);
}

MatrixXd EffFieldCalc::exchangeField(const MRef &Mag)
{
	return XC_field_OP * Mag;
}

void EffFieldCalc::setSurfaceAnisotropy()
{
	surfaceAnisotropyEnergyOffset = nodeSurfaceArea.cwiseProduct(Ksn).sum();
	Ksn = Ksn.cwiseProduct(nodeSurfaceArea);
}

double EffFieldCalc::surfaceAnisotropyEnergy(MRef &Mag)
{
	double esurf = 0.;
#ifdef _OPENMP
	double ktms;
	const int nx = Mag.rows();
#pragma omp parallel for private(ktms) reduction(+ \
												 : esurf)
	for (int i = 0; i < nx; ++i)
	{
		ktms = nv_nx.row(i).dot(Mag.row(i));
		esurf += Ksn(i) * (1. - ktms * ktms);
	}
#else
	VectorXd t1 = nv_nx.cwiseProduct(Mag).rowwise().sum();
	esurf = surfaceAnisotropyEnergyOffset - Ksn.cwiseProduct(t1).dot(t1);
#endif
	return esurf;
}

MatrixXd EffFieldCalc::surfaceAnisotropyField(MRef &Mag)
{
#ifdef _OPENMP
#pragma omp parallel shared(Mag)
	{
		int chunk_size = (nx / omp_get_num_threads());
#pragma omp for schedule(static, chunk_size)
#endif
		for (int i = 0; i < nx; ++i)
		{
			double nm = Mag.row(i).dot(nv_nx.row(i));
			Hsrf.row(i) = -2. * Ksn(i) * nm * (nm * Mag.row(i) - nv_nx.row(i)) * invNodeVol(i);
		}
#ifdef _OPENMP
	}
#endif
	// #else
	// 	  Hsrf = 2. * Ksn.cwiseProduct(invNodeVol).cwiseProduct(nv_nx.cwiseProduct(Mag).rowwise().sum()).asDiagonal() * nv_nx ;
	// #endif
	// //	  Hsrf += surfIntKs(Mag); // unnecessary, since RW field at boundary is collinear with m and thus does not affect the LLG dynamics
	return Hsrf;
}

/*
MatrixXd EffFieldCalc::surfIntKs( MRef &Mag ) {
// Implementing Rado-Weertman BCs
// Note that Ksn is already multiplied with NodeSurfaceArea and invNodeVol
	return  Ksn.asDiagonal() * Mag;
}
*/

MatrixXd EffFieldCalc::surfIntDMI(MRef &Mag)
{
	// Boundary conditions according to Rohart and Thiaville
	MatrixXd surfDMI = nodeSurfaceArea.asDiagonal() * cross(Mag, nv_nx);
	return surfDMI;
}

MatrixXd EffFieldCalc::dmiField(MRef &Mag)
{
	Hdmi = (surfIntDMI(Mag) - 2. * curlM(Mag)).array().colwise() * invNodeVol.cwiseProduct(D).array();
	return Hdmi;
}

double EffFieldCalc::dmiEnergy(MRef &Mag)
{
	return (NodeVol.transpose() * Hdmi.cwiseProduct(Mag)).sum() / (-2.);
}

double EffFieldCalc::exchEnergy_direct(MRef &Mag){
	VectorXd mx = Mag.col(x);
	VectorXd my = Mag.col(y);
	VectorXd mz = Mag.col(z);
	int nx = mx.size();
	VectorXd dmxdx = gradX * mx;
	VectorXd dmxdy = gradY * mx;
	VectorXd dmxdz = gradZ * mx;
	VectorXd dmydx = gradX * my;	
	VectorXd dmydy = gradY * my;
	VectorXd dmydz = gradZ * my;
	VectorXd dmzdx = gradX * mz;
	VectorXd dmzdy = gradY * mz;	
	VectorXd dmzdz = gradZ * mz;

	VectorXd exch(nx);
	for (int i = 0; i < nx; ++i) {
		exch(i) = dmxdx(i)*dmxdx(i) + dmxdy(i)*dmxdy(i) + dmxdz(i)*dmxdz(i) +
  				  dmydx(i)*dmydx(i) + dmydy(i)*dmydy(i) + dmydz(i)*dmydz(i) + 
				  dmzdx(i)*dmzdx(i) + dmzdy(i)*dmzdy(i) + dmzdz(i)*dmzdz(i);
		exch(i) *= A(i) * NodeVol(i);
	}
	return exch.sum() / NodeVol.sum();	


}

double EffFieldCalc::dmiEnergy_direct(MRef &Mag)
{
	VectorXd mx = Mag.col(x);
	VectorXd my = Mag.col(y);
	VectorXd mz = Mag.col(z);
	int nx = mx.size();

	VectorXd dmxdx = gradX * mx;
	VectorXd dmxdy = gradY * mx;
	VectorXd dmxdz = gradZ * mx;
	VectorXd dmydx = gradX * my;	
	VectorXd dmydy = gradY * my;
	VectorXd dmydz = gradZ * my;
	VectorXd dmzdx = gradX * mz;
	VectorXd dmzdy = gradY * mz;
	VectorXd dmzdz = gradZ * mz;

	MatrixXd curlM(nx,3);
	curlM.col(x) = dmzdy - dmydz;
	curlM.col(y) = dmxdz - dmzdx;
	curlM.col(z) = dmydx - dmxdy;

	VectorXd mCurlm = VectorXd::Zero(nx);
	for (int i = 0; i < nx; ++i) {
		mCurlm(i) = Mag.row(i).dot(curlM.row(i));
	}
	return mCurlm.transpose() * D.cwiseProduct(NodeVol);
}


double EffFieldCalc::cubicAnisotropyEnergy(MRef &Mag)
{
	double Ecub = 0.;
	Vector3d alf;
#ifdef _OPENMP
#pragma omp parallel for private(alf) shared(Mag) reduction(+ \
															: Ecub)
#endif
	for (int i = 0; i < nx; ++i)
	{
		alf = cubicAxes[i].transpose() * Mag.row(i).transpose(); // note: easy axes are stored as columns in cubicAxes[i]
		alf = alf.cwiseProduct(alf);
		Ecub += NodeVol(i) * (Kc1(i) * (alf(0) * alf(1) + alf(1) * alf(2) + alf(0) * alf(2)) + Kc2(i) * alf.prod());
	}
	return Ecub;
}

MatrixXd EffFieldCalc::cubicAnisotropyField(MRef &Mag)
{
	Vector3d dea, alf;
	double a0, a1, a2;
#ifdef _OPENMP
#pragma omp parallel for private(alf, dea, a0, a1, a2) shared(Mag)
#endif
	for (int i = 0; i < nx; ++i)
	{
		alf = cubicAxes[i].transpose() * Mag.row(i).transpose(); // this is a 3x3 matrix-vector product
		// use first cubicAxes[i].transpose() because easy axes are stored in cols, second one to convert Mag.row(i) into column vector for MVP
		for (int j = 0; j < 3; ++j)
		{
			a0 = alf(j);
			a1 = alf((j + 1) % 3);
			a2 = alf((j + 2) % 3);
			dea(j) = -2. * a0 * (Kc1(i) * (a1 * a1 + a2 * a2) + Kc2(i) * a1 * a1 * a2 * a2);
		}
		Hcub.row(i) = cubicAxes[i] * dea;
		// here cubicAxes[i] is the Jacobi determinant [du/dm]
		// where u1, u2, u3 are the magnetization components along the cubic axes
	}
	return Hcub;
}

void EffFieldCalc::setDMIField(const MRef &Hdmi_)
{
	// this is needed to transfer the GPU-calculated field so that the DMI energy can be calculated on the host.
	Hdmi = Hdmi_;
}

MatrixXd EffFieldCalc::curlM(MRef &Mag)
{
	MatrixXd curlM(nx, 3);
	/*
// copying to const Vector doesn't affect the performance.
	const VectorXd Mx = Mag.col(x);
	const VectorXd My = Mag.col(y);
	const VectorXd Mz = Mag.col(z);
	curlM.col(x) = tGradY * Mz - tGradZ * My;
	curlM.col(y) = tGradZ * Mx - tGradX * Mz;
	curlM.col(z) = tGradX * My - tGradY * Mx;
*/

	curlM.col(x) = tGradY * Mag.col(z) - tGradZ * Mag.col(y);
	curlM.col(y) = tGradZ * Mag.col(x) - tGradX * Mag.col(z);
	curlM.col(z) = tGradX * Mag.col(y) - tGradY * Mag.col(x);

	return curlM;
}
