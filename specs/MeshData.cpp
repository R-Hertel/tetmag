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
 * MeshData.cpp
 *
 *  Created on: May 2, 2017
 *      Author: riccardo
 */

#include "MeshData.h"
#include "BEMprocessing.h"
#include <iostream>
#include <Eigen/Sparse>
#include <string>
#include <ios>
#include <fstream>
#include "Materials.h"
#include "auxiliaries.h"

using namespace Eigen;
enum coords {x, y, z};
enum mats {one, two, three};

void MeshData::ascribeMaterialsToNodes() {
//	if (fileExists("nodes.mat")) {std::cout << "reading material assignment from file 'nodes.mat'" << std::endl;}
//	exit(1);
	Material mat;
//	int nx = xyz.rows();
//	int numberOfMaterials;
//	NodeMaterial = VectorXi::Zero(nx);
	NodeMaterial = mat.assignMaterials(xyz);
/*
 * Ks selection for Co wedge; only on bottom surface and with different values on each half of the surface
 *
	for (int i = 0; i < nx; ++i) {
		NodeMaterial(i) = one; // positive Ks
		if (xyz(i, y) > 0.25) NodeMaterial(i) = three; // negative Ks on opposite domain along x
		if (nv_nx(i, z) > -0.9) {
			NodeMaterial(i) = two; // set Ks to zero if not on bottom surface
		}
	}
*/
	// If needed, insert position-dependent assignment of material type here
}


void MeshData::selectAnisotropicSurfaces(const VectorXd& Ks) {

(void) Ks;	// only to silence the "unused parameter" warning when nothing is selected

//	Deselect the surfaces where surface anisotropy does not apply.
//  The selection can be done using the node position xyz and/or the direction of the node normal vector nv_nx.
//  Set Ks(i) to zero in the corresponding nodes i.
//  Note: Node numbering is global, modified Ks vector will be returned by reference.
//
//  Example 1: select only surfaces located at z >= 0.4
//	int nx = xyz.rows();
//	for (int i = 0; i < nx; ++i) {
//		if (xyz(i,z) < 0.4) Ks(i) = 0.0;
//	}

//	Example 2: Select only surfaces with y-component of the normal vector > 0.8
//	int nx = xyz.rows();
//	for (int i = 0; i < nx; ++i) {
//		if (nv_nx(i, y) <= 0.8)
//			Ks(i) = 0.0;
//	}
}


void MeshData::preprocessBEM(bool useH2) {
	BEMprocessing bem(*this);
	bem.preProcess(useH2, name);
	boundaryNodes = bem.getBoundaryNodes();
	laplaceBEM = bem.getLaplaceBEMmatrix();
	nv_nx = bem.getBoundaryNodeNormalVectors();
	nodeArea = bem.getNodeArea();
}
