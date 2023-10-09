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
	Hl = prog.getHloc();
	Hl.frequency /= 1.e3; // convert input data form GHz to 1/ps
	Hl.omega = 2. * PhysicalConstants::pi * Hl.frequency;
	Hl.amplitude /= 1000.; // convert input data from mT to T
	Hl.amplitude /= PhysicalConstants::mu0; // convert [T] to [A/m]
	useGPU = (prog.solverType == "gpu" || prog.solverType == "pl");
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


void SimulationData::readLocalFieldProfile() {
	if (!Hl.isUsed) return;
	std::string profileFileName = "fieldProfile.vtu";
	if (!fileExists(profileFileName)) {
		std::cout << "ERROR: Option for local field excitation selected, but no file with local field profile found.\n";
		std::cout << "Prepare a file named \"fieldProfile.vtu\" with the field profile information, or deselect the option." << std::endl;
		exit(1);
		return;
	}
	VTKReader r(profileFileName);
	MatrixXd Profile = r.readMag();
	bool goodNodes = r.checkNodeNumber(nx);
	if (!goodNodes) {
		std::cout << "Please provide a compatible '" << profileFileName << "' file, or deselect the the 'local field' option." << std::endl;
		exit(1);
	}

	// ToDo: replace this with a more user-friendly file reader / input routine for the "profile" option. //
	double maxVal = Profile.rowwise().norm().maxCoeff();
	if (areEqual(maxVal, 0.))  {
		std::cout << "ZERO field profile detected. Not using local field." << std::endl;
	    Hl.isUsed = false;
	    return;
	}
	Profile /= maxVal;
	Hl.setProfile(Profile);
	std::cout << "Using local field profile." << std::endl;
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
