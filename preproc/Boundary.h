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
 * MeshBoundary.h
 *
 *  Created on: Oct 13, 2016
 *      Author: riccardo
 */

#ifndef BOUNDARY_H_
#define BOUNDARY_H_

#include <Eigen/Dense>
#include <vector>
#include "BEProperties.h"

class MeshBoundary {
private:
	Eigen::MatrixXd xyz;
	Eigen::MatrixXi fel;
	Eigen::VectorXi boundaryNodeIndex;
	unsigned int nbel;
	unsigned int bnx;
	std::vector<TriangleFace> BoundaryTriangle;
	std::vector<int> boundaryNodes;
	double totalSurfaceByNodes;
	double totalSurfaceByElements;
	void findBoundaryElements();
	Eigen::VectorXd ElementSurface;
	void calculateElementSurface();
	void findBoundaryNodes();
	Eigen::VectorXd nodeSurfaceArea;
	void calculateBoundaryNodeSurface();
	std::vector<BEproperties> belprops;
	Eigen::MatrixX3d nv_el;
	Eigen::MatrixX3d nv_n;
	void calculateNormalVector();
	Eigen::Vector3d CenterOfMass(const Eigen::MatrixXi&, const Eigen::MatrixXd&, int );
	void prepareBoundaryNodeIndex();
	void assignNormalVectorsToNodes();
public:
	MeshBoundary(const Eigen::MatrixXd&, const Eigen::MatrixXi&);
	void extractFeatures();
	std::vector<TriangleFace> getBoundaryElements();
	Eigen::VectorXd getElementSurface();
	double getTotalSurfaceByElements();
	std::vector<int> getBoundaryNodes();
	double getTotalSurfaceByNodes();
	Eigen::MatrixX3d getNormalVectors();
	std::vector<BEproperties> getBELproperties();
	Eigen::VectorXi getBoundaryNodeIndex();
	Eigen::MatrixX3d getNodeNormalVectors();
	Eigen::VectorXd getNodeArea();
};

#endif /* BOUNDARY_H_ */
