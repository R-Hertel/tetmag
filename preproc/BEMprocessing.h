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
 * BEMProcessing.h
 *
 *  Created on: Oct 11, 2016
 *      Author: riccardo
 */

#ifndef BEMPROCESSING_H_
#define BEMPROCESSING_H_

#include "MeshData.h"
#include <Eigen/Dense>
#include <vector>
#include <string>
#include "BEProperties.h"
#include "Lindholm.h"

struct tet_edge;
class BEMprocessing {
private:
	MeshData msh;
	size_t nx;
	size_t nel;
	size_t bnx;
	size_t nbel;
	Eigen::VectorXd sub;
	Eigen::VectorXi boundaryNodeIndex; // maps global node number to boundary node number
	Eigen::MatrixXi el_ed;
	Eigen::MatrixXi bel_b;
	Eigen::MatrixXd bxyz;
	Eigen::MatrixXd nv_n;
	Eigen::VectorXd nodeArea;
	void findBoundaryElementsAndNodes();
	void calculateSolidAngles();
	void assembleLaplaceBEMmatrix();
	Eigen::MatrixXd laplaceBEMmat;
	std::vector<BEproperties> belprops;
	void testLaplaceMatrix();
	bool matrixExists();
	void storeLaplaceBEMMatrix();
	Eigen::MatrixXd readLaplaceBEMMatrix();
	void HLibSuite();
	void H2LibCSuite(std::string);
	std::vector<tet_edge> prepareEdges();
	void reorderNodesCCW();
	void outputBoundary();
	void prepareLocalBoundaryMatrices();
	Lindholm lindholm;
public:
	BEMprocessing(const MeshData&);
	void preProcess(bool, std::string);
	std::vector<int> getBoundaryNodes();
	Eigen::MatrixXd getLaplaceBEMmatrix();
	Eigen::MatrixXd getBoundaryNodeNormalVectors();
	Eigen::VectorXd getNodeArea();
};

#endif /* BEMPROCESSING_H_ */
