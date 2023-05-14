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
 * MeshData.h
 *
 *  Created on: May 2, 2017
 *      Author: riccardo
 */

#ifndef MESHDATA_H_
#define MESHDATA_H_
#include <Eigen/Dense>
#include <vector>
#include <string>
#include <Eigen/SparseCore>
#include "typedefs.h"


class MeshData {
public:
	std::string name;
	Eigen::MatrixXi fel;
	Eigen::MatrixXd xyz;
	Eigen::VectorXd NodeVolume;
	std::vector<int> boundaryNodes;
	Eigen::MatrixXd nv_nx; // normalized surface anisotropy axis vector, Defined at each node. Zero for non-boundary nodes.
	Eigen::VectorXd nodeArea; // area assigned to boundary nodes. Defined for all nodes. Zero for non-boundary nodes.
	SpMat stiff;
	SpMat gradX;
	SpMat gradY;
	SpMat gradZ;
	SpMat tGradX;
	SpMat tGradY;
	SpMat tGradZ;
	SpMat dirichletMatrix;
	Eigen::VectorXi NodeMaterial;
	void ascribeMaterialsToNodes();
	void selectAnisotropicSurfaces(const Eigen::VectorXd&);
	Eigen::MatrixXd laplaceBEM; // dense BEM matrix, only used in the case prog.useH2 = false
	void preprocessBEM(bool);
};

#endif /* MESHDATA_H_ */
