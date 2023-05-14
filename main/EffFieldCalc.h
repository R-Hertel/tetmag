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
 * EffFieldCalc.h
 *
 *  Created on: May 25, 2021
 *      Author: riccardo
 */

#ifndef MAIN_EFFFIELDCALC_H_
#define MAIN_EFFFIELDCALC_H_
#include <Eigen/SparseCore>
#include <Eigen/Core>
#include "typedefs.h"
#include "SimulationData.h"
#include "MeshData.h"
#include <vector>

class EffFieldCalc {
	int nx;
	SpMat XC_field_OP;
	Eigen::VectorXd A;
	Eigen::VectorXd D;
	Eigen::VectorXd Ku1;
	Eigen::MatrixXd KuAxis;
	std::vector<Eigen::Matrix3d> cubicAxes;
	Eigen::VectorXd Kc1, Kc2;
	Eigen::VectorXd Ksn; // Ks at nodes, zero for non-boundary nodes.
	Eigen::VectorXd nodeSurfaceArea;
	Eigen::MatrixXd nv_nx; // normal surface vector, defined at all nodes. Zero vector for non-boundary nodes. Normalized.
	SpMat tGradX;
	SpMat tGradY;
	SpMat tGradZ;
	SpMat gradX;
	SpMat gradY;
	SpMat gradZ;	
	Eigen::MatrixXd Hexc;
	Eigen::MatrixXd Hani;
	Eigen::MatrixXd Hsrf;
	Eigen::MatrixXd Hdmi;
	Eigen::MatrixXd Hcub;
	Eigen::VectorXd NodeVol;
	Eigen::VectorXd invNodeVol;

	double uniaxialAnisotropyEnergyOffset;
	double surfaceAnisotropyEnergyOffset;

	void setExchangeFieldOperator();
	void setUniaxialAnisotropy();
	void setSurfaceAnisotropy();
	Eigen::MatrixXd curlM(MRef&);
	Eigen::MatrixXd surfIntDMI(MRef&);
public:
	EffFieldCalc(SimulationData&, const MeshData&);

	Eigen::MatrixXd exchangeField(const MRef&);
	double calcExchangeEnergy(const MRef&);
	double exchEnergy_direct(MRef&);
	Eigen::MatrixXd uniaxialAnisotropyField(const MRef&);
	double uniaxialEnergy(const MRef&);
	Eigen::MatrixXd surfaceAnisotropyField(MRef&);
	double surfaceAnisotropyEnergy(MRef&);
	Eigen::MatrixXd dmiField(MRef&);
	double dmiEnergy(MRef&);
	double dmiEnergy_direct(MRef&);
	Eigen::MatrixXd cubicAnisotropyField(MRef&);
	double cubicAnisotropyEnergy(MRef&);
	void setDMIField(const MRef&);
};

#endif /* MAIN_EFFFIELDCALC_H_ */
