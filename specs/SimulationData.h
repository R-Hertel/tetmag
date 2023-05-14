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
 * SimulationData.h
 *
 *  Created on: Apr 27, 2017
 *      Author: riccardo
 */

#ifndef SIMULATIONDATA_H_
#define SIMULATIONDATA_H_
#include <string>
#include <vector>
#include <Eigen/Dense>
#include "ProgramSpecs.h"
#include "MeshData.h"
#include "Materials.h"

class SimulationData {
	void allocateMaterialVectors();
	void prepareUniaxialAnisotropyAtNodes(MeshData& , std::vector<Material>& );
public:
	int nx;
	int bnx;
	Eigen::VectorXd Ks;   // Surface anisotropy defined at each node. Zero for non-surface nodes.
	Eigen::MatrixXd Kuni; // vector of uniaxial easy axis
	Eigen::VectorXd Ku1;  // first uniaxial anisotropy constant
	Eigen::VectorXd Js;
	Eigen::VectorXd Kc1, Kc2;
	std::vector<Eigen::Matrix3d> cubicAxes;
	Eigen::VectorXd A;
	Eigen::VectorXd D;    // DMI constant (bulk)
	double scale;
	double alpha;
	double Hext;
	double theta_H;
	double phi_H;
	bool Hys = false;
	double H_Hys = 0., theta_Hys = 0., phi_Hys = 0.;
	double psTimeScaleFactor();
	void getProgramData(ProgramSpecs&);
	void scaleToRealSize();
	void setMaterialParametersAtNodes(MeshData&, std::vector<Material>& );
	Hdynamic Hp;
	Hlocal Hl;
	void readLocalFieldProfile();
	bool useGPU;
};

#endif /* SIMULATIONDATA_H_ */
