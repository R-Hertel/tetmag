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
 * PreprocessMesh.h
 *
 *  Created on: Oct 18, 2016
 *      Author: riccardo
 */

#ifndef FEMPROCESSING_H_
#define FEMPROCESSING_H_

#include "MeshData.h"
#include <Eigen/SparseCore>
#include <Eigen/Dense>
#include <list>
#include <vector>
#include <cstddef>
#include <string>
#include "typedefs.h"

class PreprocessMesh {
private:
	size_t nx;
	size_t ntet;
//	typedef Eigen::SparseMatrix<double, Eigen::RowMajor> SpMat;
	Eigen::MatrixXd ShapeFuncGradX;
	Eigen::MatrixXd ShapeFuncGradY;
	Eigen::MatrixXd ShapeFuncGradZ;
	SpMat gradOpX;
	SpMat gradOpY;
	SpMat gradOpZ;
	SpMat tGradOpX;
	SpMat tGradOpY;
	SpMat tGradOpZ;
	Eigen::VectorXd ElementVolume;
	SpMat stiff;
	std::vector<int> DirichletNodes;
	std::vector<double> DirichletValues;
	void calculateVolumes();
	void calcShapeFuncGrads();
	void calculateTransposeGradientOperators();
	void calculateGradientOperators();
	double signedTetrahedronVolume6(int);
	double signedTetrahedronVolume6(const Eigen::Vector3d&, const Eigen::Vector3d&, const Eigen::Vector3d&, const Eigen::Vector3d& );
	void testGradientOperators();
	Eigen::Matrix4d elementStiffnessMatrix(int);
	Eigen::Matrix4d elementMassMatrix(int);
	void assembleStiffnessMatrix();
	void addElementMatrixToGlobalMatrix(SpMat&, const Eigen::Matrix4d&, int);
	MeshData msh;
	void prepareDirichletMatrix();
public:
	PreprocessMesh(MeshData&);
	void prepareAllMatrices();
	MeshData getMeshData();
	Eigen::VectorXi getDirichletNodes();
	Eigen::VectorXd getDirichletValues();
	void checkVolumes();
};

#endif /* FEMPROCESSING_H_ */
