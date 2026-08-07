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
 * SimulationData.cpp
 *
 *  Created on: Apr 27, 2017
 *      Author: riccardo
 */

#include "SimulationData.h"
#include "ProgramSpecs.h"
#include "PhysicalConstants.h"
#include "MeshData.h"
#include <vector>
#include "Materials.h"
#include <Eigen/Dense>
#include <iostream>
#include "auxiliaries.h"
#include "VTKReader.h"

using namespace Eigen;

void SimulationData::scaleToRealSize(){
	A /= (scale * scale);
	Ks /= scale;
	D /= scale;
}


double SimulationData::psTimeScaleFactor() {
	const double pico = 1.e-12;
//	double scf = PhysicalConstants::gamma0 / (1. + alpha * alpha) * pico; // gamma in [ m / (A.s) ]
//	double scf = PhysicalConstants::gamma0 * PhysicalConstants::mu0 / (1. + alpha * alpha) * pico; // gamma in [ 1 / (T.s) ]
	return gamma / (1. + alpha * alpha) * pico;
}


void SimulationData::getProgramData(ProgramSpecs& prog) {
	alpha = prog.alpha;
	Hext = prog.getExternalField();
	Hext /= 1000.; // convert input data from mT to T
	theta_H = std::get<0> ( prog.getHextAngles() );
	phi_H   = std::get<1> ( prog.getHextAngles() );
	scale = prog.scale;
	Hp = prog.getHdyn();
	Hp.pulsePeak /= 1000.; // convert input data from mT to T
	Hp.sweepStart /= 1000.;
	Hp.sweepEnd /= 1000.;
	Hp.rfAmplitude /= 1000.; // convert input data from mT to T
	Hp.rfFrequency /= 1.e3;  // convert input data from GHz to 1/ps
	Hp.rfOmega = 2. * PhysicalConstants::pi * Hp.rfFrequency;
	Hp.staticLocalAmplitude /= PhysicalConstants::mu0; // convert T -> A/m
	useGPU = (prog.solverType == "gpu");
	gamma = prog.gamma;

}


void SimulationData::allocateMaterialVectors() {
	VectorXd zeroVec = VectorXd::Zero(nx);
	Ku1 = zeroVec;
	Js = zeroVec;
	Kc1 = zeroVec;
	Kc2 = zeroVec;
	cubicAxes.resize(nx);
	Ks = zeroVec;
	A = zeroVec;
	D = zeroVec;
}


void SimulationData::readFieldProfile() {
	bool needProfile = Hp.pulseHasProfile || Hp.rfHasProfile || Hp.staticLocalIsUsed;
	if (!needProfile) return;
	std::string profileFileName = "fieldProfile.vtu";
	if (!fileExists(profileFileName)) {
		std::cerr << "ERROR: Spatial profile option selected, but 'fieldProfile.vtu' not found." << std::endl;
		std::cerr << "Provide the file or disable profile options." << std::endl;
		exit(1);
	}
	VTKReader r(profileFileName);
	MatrixXd Profile = r.readMag();
	bool goodNodes = r.checkNodeNumber(nx);
	if (!goodNodes) {
		std::cerr << "Node count mismatch in 'fieldProfile.vtu'." << std::endl;
		exit(1);
	}
	double maxVal = Profile.rowwise().norm().maxCoeff();
	if (areEqual(maxVal, 0.)) {
		std::cerr << "ZERO field profile detected. Disabling profile options." << std::endl;
		Hp.pulseHasProfile   = false;
		Hp.rfHasProfile      = false;
		Hp.staticLocalIsUsed = false;
		return;
	}
	Profile /= maxVal;
	Hp.fieldProfile = Profile;
	if (Hp.staticLocalIsUsed)
		Hp.staticLocalH = Profile * Hp.staticLocalAmplitude;
	std::cout << "Using spatial field profile from '" << profileFileName << "'." << std::endl;
}


void SimulationData::setMaterialParametersAtNodes(MeshData& msh, std::vector<Material>& mats) {
	nx = msh.xyz.rows();
	allocateMaterialVectors();
	prepareUniaxialAnisotropyAtNodes(msh, mats);
	int i, j;
#ifdef _OPENMP
#pragma omp parallel for private(i, j)
#endif
	for (i = 0; i < nx; i++) {
		j = msh.NodeMaterial(i);
		Js(i)  = mats[j].Js;
		cubicAxes[i] = mats[j].calculateCubicAxesVectors();
		Ku1(i) = mats[j].Ku  ;
		Kc1(i) = mats[j].Kc1 ;
		Kc2(i) = mats[j].Kc2 ;
		Ks(i)  = mats[j].Ks  ;
		A(i)   = mats[j].A   ;
		D(i)   = mats[j].D   ;
	}
}


void SimulationData::prepareUniaxialAnisotropyAtNodes(MeshData& msh, std::vector<Material>& mats) {
	enum coords {x, y, z};
	double gradToRad = PhysicalConstants::pi / 180.;
	Kuni = MatrixXd::Zero(nx,3);
	int i, j;
#ifdef _OPENMP
#pragma omp parallel for private(i, j)
#endif
	for (i = 0; i < nx; i++) {
		j = msh.NodeMaterial(i);
		double theta = mats[j].theta_u * gradToRad;
		double phi = mats[j].phi_u * gradToRad;
		Kuni(i, x) = std::sin(theta) * std::cos(phi);
		Kuni(i, y) = std::sin(theta) * std::sin(phi);
		Kuni(i, z) = std::cos(theta);
	}
}
