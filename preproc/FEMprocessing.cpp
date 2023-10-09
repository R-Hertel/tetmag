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
 * PreprocessMesh.cpp
 *
 *  Created on: Oct 18, 2016
 *      Author: riccardo
 */

#include "FEMprocessing.h"
#include "auxiliaries.h"
#include <iostream>
#include <fstream>
#include <cstddef> 		// for size_t
#include <cmath> 		// for pow()
#include <Eigen/Dense>
#include <Eigen/SparseLU>
#include <cassert>
#include <vector>

#include <string>
#include <sstream>
#include <ios>
#include <fstream>
#include "Boundary.h"

#include <cstdio>
#include <memory>
#include <stdexcept>
#include <array>

using namespace Eigen;


PreprocessMesh::PreprocessMesh(MeshData& msh_) :
		nx(msh_.xyz.rows()),
		ntet(msh_.fel.rows()),
		gradOpX(nx, nx),// <-- This needs to be in the constructor, not in the header file.
		gradOpY(nx, nx),
		gradOpZ(nx, nx),
		tGradOpX(nx, nx),
		tGradOpY(nx, nx),
		tGradOpZ(nx, nx),
		stiff(nx, nx),
//		mass(nx, nx),
		msh(msh_){}

void PreprocessMesh::prepareAllMatrices() {
	calculateVolumes();
	calcShapeFuncGrads();
	calculateGradientOperators();
	gradOpX.makeCompressed();
	gradOpY.makeCompressed();
	gradOpZ.makeCompressed();

	tGradOpX.makeCompressed();
	tGradOpY.makeCompressed();
	tGradOpZ.makeCompressed();

	assembleStiffnessMatrix();
	stiff.makeCompressed();

	msh.gradX = gradOpX;
	msh.gradY = gradOpY;
	msh.gradZ = gradOpZ;
	msh.tGradX = tGradOpX;
	msh.tGradY = tGradOpY;
	msh.tGradZ = tGradOpZ;
	msh.stiff = stiff;
	prepareDirichletMatrix();
}


void PreprocessMesh::checkVolumes() {
	std::cout << "Total volume: " << ElementVolume.sum() << std::endl;
	std::cout << "average volume: " << ElementVolume.sum() / ntet << std::endl;
	std::cout << "Smallest volume: " << ElementVolume.minCoeff() << std::endl;
	std::cout << "Largest volume: " << ElementVolume.maxCoeff() << std::endl;
}


void PreprocessMesh::calculateVolumes() {
	int nel = msh.fel.rows();
	ElementVolume  = VectorXd::Zero(nel);
	msh.NodeVolume = VectorXd::Zero(nx);
	for (int i = 0; i < nel; ++i) {
		ElementVolume(i) = std::abs(signedTetrahedronVolume6(i)) / 6.;
		for (int j = 0; j < 4; ++j) {
			msh.NodeVolume(msh.fel(i, j)) += ElementVolume(i) / 4.;}
	}
}

double PreprocessMesh::signedTetrahedronVolume6(int el) {
	// Tetrahedron volume multiplied by 6, including the sign of the determinant.
	Matrix4d el_mat;
	for (int n = 0; n < 4; ++n) {
		int global_n = msh.fel(el, n);
		el_mat.row(n) << 1.0, msh.xyz.row(global_n);
	}
	return el_mat.determinant();
}


MeshData PreprocessMesh::getMeshData(){ return msh; }


VectorXd PreprocessMesh::getDirichletValues() {
	VectorXd EigenDirichletValues = Map<VectorXd>(DirichletValues.data(), DirichletValues.size());
	return EigenDirichletValues;
}


VectorXi PreprocessMesh::getDirichletNodes(){
	VectorXi EigenDirichletNodes = Map<VectorXi>(DirichletNodes.data(), DirichletNodes.size());
	return EigenDirichletNodes;
}


double PreprocessMesh::signedTetrahedronVolume6(const Vector3d& p1, const Vector3d& p2, const Vector3d& p3, const Vector3d& p4) {
	Matrix4d m;
	enum coords {x, y, z};
	m   <<	  1.,	 1., 	1.,    1.,
			p1(x), p2(x), p3(x), p4(x),
			p1(y), p2(y), p3(y), p4(y),
			p1(z), p2(z), p3(z), p4(z);
	double V6 = m.determinant();
	return V6;
}


void PreprocessMesh::calcShapeFuncGrads() {
	// cf. Zienkewicz et al., "The Finite Element Method", 6th ed.,
	// section 4.11, "Tetrahedral Elements", p. 123.
	const int dirs = 3;
	ShapeFuncGradX = MatrixXd::Zero(ntet, 4);
	ShapeFuncGradY = MatrixXd::Zero(ntet, 4);
	ShapeFuncGradZ = MatrixXd::Zero(ntet, 4);
	MatrixXd Element(4,3);
	Matrix3d node_mat;
	enum Coords { x, y, z };
	enum Points { node1, node2, node3, node4 };
	size_t el;
#ifdef _OPENMP
#pragma omp parallel for
#endif
	for (el = 0; el < ntet; ++el) {
#ifdef _OPENMP
#pragma	omp critical
#endif
		{
		for (int i = 0; i < 4; i++) {
			Element.row(i) = msh.xyz.row(msh.fel(el, i));
		}
		double V6 = signedTetrahedronVolume6(
				Element.row(node1), Element.row(node2),	Element.row(node3), Element.row(node4));
		for (int n = 0; n < 4; ++n) {
			for (int i = 0; i < dirs; ++i) {
				node_mat.row(i) = Element.row((n + 1 + i) % 4); // cyclic permutation
			}
			double grad_sign = pow(-1., n + 1); 
			Matrix3d gradMatX = node_mat;
			Matrix3d gradMatY = node_mat;
			Matrix3d gradMatZ = node_mat;
			gradMatX.col(x).setOnes();
			gradMatY.col(y).setOnes();
			gradMatZ.col(z).setOnes();
			ShapeFuncGradX(el, n) = grad_sign * gradMatX.determinant() / V6;
			ShapeFuncGradY(el, n) = grad_sign * gradMatY.determinant() / V6;
			ShapeFuncGradZ(el, n) = grad_sign * gradMatZ.determinant() / V6;
			}
		}
	}
}

void PreprocessMesh::calculateTransposeGradientOperators() {
// transfered to calcGradOp
/*
	double node_weight;
	gradOpX.reserve(VectorXi::Constant(nx, 30));
	gradOpY.reserve(VectorXi::Constant(nx, 30));
	gradOpZ.reserve(VectorXi::Constant(nx, 30));

	tGradOpX.reserve(VectorXi::Constant(nx, 30));
	tGradOpY.reserve(VectorXi::Constant(nx, 30));
	tGradOpZ.reserve(VectorXi::Constant(nx, 30));
	for (size_t el = 0; el < ntet; el++) {
		for (int local_i = 0; local_i < 4; ++local_i) {
			for (int local_j = 0; local_j < 4; ++local_j) {
				int global_i = msh.fel(el, local_i);
				int global_j = msh.fel(el, local_j);
				node_weight = ElementVolume(el) / 4.; // <- Volume integration over one of the shape functions
				tGradOpX.coeffRef(global_i, global_j) += ShapeFunctionGradientX(el, local_i) * node_weight;
				tGradOpY.coeffRef(global_i, global_j) += ShapeFunctionGradientY(el, local_i) * node_weight;
				tGradOpZ.coeffRef(global_i, global_j) += ShapeFunctionGradientZ(el, local_i) * node_weight;

				gradOpX.coeffRef(global_i, global_j) += ShapeFunctionGradientX(el, local_j) * node_weight / msh.NodeVolume(global_i);
				gradOpY.coeffRef(global_i, global_j) += ShapeFunctionGradientY(el, local_j) * node_weight / msh.NodeVolume(global_i);
				gradOpZ.coeffRef(global_i, global_j) += ShapeFunctionGradientZ(el, local_j) * node_weight / msh.NodeVolume(global_i);
			}
		}
	}
*/
}


void PreprocessMesh::calculateGradientOperators() {
	double node_weight;
	gradOpX.reserve(VectorXi::Constant(nx, 30));
	gradOpY.reserve(VectorXi::Constant(nx, 30));
	gradOpZ.reserve(VectorXi::Constant(nx, 30));

	tGradOpX.reserve(VectorXi::Constant(nx, 30));
	tGradOpY.reserve(VectorXi::Constant(nx, 30));
	tGradOpZ.reserve(VectorXi::Constant(nx, 30));
	for (size_t el = 0; el < ntet; el++) {
		for (int local_i = 0; local_i < 4; ++local_i) {
			for (int local_j = 0; local_j < 4; ++local_j) {
				int global_i = msh.fel(el, local_i);
				int global_j = msh.fel(el, local_j);
				node_weight = ElementVolume(el) / 4.; // <- Volume integration over one of the shape functions
				gradOpX.coeffRef(global_i, global_j) += ShapeFuncGradX(el, local_j) * node_weight / msh.NodeVolume(global_i);
				gradOpY.coeffRef(global_i, global_j) += ShapeFuncGradY(el, local_j) * node_weight / msh.NodeVolume(global_i);
				gradOpZ.coeffRef(global_i, global_j) += ShapeFuncGradZ(el, local_j) * node_weight / msh.NodeVolume(global_i);

				tGradOpX.coeffRef(global_i, global_j) += ShapeFuncGradX(el, local_i) * node_weight;
				tGradOpY.coeffRef(global_i, global_j) += ShapeFuncGradY(el, local_i) * node_weight;
				tGradOpZ.coeffRef(global_i, global_j) += ShapeFuncGradZ(el, local_i) * node_weight;
			}
		}
	}
}


void PreprocessMesh::addElementMatrixToGlobalMatrix(SpMat& globalMat, const Matrix4d& elementMat, int el_num) {
	for (int local_i = 0; local_i < 4; ++local_i) {
		for(int local_j = 0; local_j < 4; ++local_j) {
			int global_i = msh.fel(el_num, local_i);
			int global_j = msh.fel(el_num, local_j);
			globalMat.coeffRef(global_i, global_j) += elementMat(local_i, local_j);
		}
	}
}

Matrix4d PreprocessMesh::elementStiffnessMatrix(int el) {
	Matrix4d elStiffMat;
	double Vel = ElementVolume(el);
	for (int i = 0; i < 4; ++i) {
		for (int j = 0; j < 4; ++j) {
			elStiffMat(i, j) = ShapeFuncGradX(el, i) * ShapeFuncGradX(el, j)
							 + ShapeFuncGradY(el, i) * ShapeFuncGradY(el, j)
							 + ShapeFuncGradZ(el, i) * ShapeFuncGradZ(el, j);
		}
	}
	elStiffMat *= Vel;
	return elStiffMat;
}


void PreprocessMesh::assembleStiffnessMatrix() {
	stiff.reserve(VectorXi::Constant(nx, 30));
	for (size_t el = 0; el < ntet; ++el) {
		addElementMatrixToGlobalMatrix(stiff, elementStiffnessMatrix(el), el);
	}
}

void PreprocessMesh::testGradientOperators() {
	VectorXd testfun(nx);
	enum coords	: int { x, y, z };
	double test_x = 0.5;
	double test_y = -17.;
	double test_z =  13.;
	double offset =  42.;
	for (size_t i = 0; i < nx; ++i) {
		double x_pos = msh.xyz(i, x);
		double y_pos = msh.xyz(i, y);
		double z_pos = msh.xyz(i, z);
		testfun(i) = test_x * x_pos + test_y * y_pos + test_z * z_pos + offset;
	}
	VectorXd fun_grad_x = gradOpX * testfun;
	VectorXd fun_grad_y = gradOpY * testfun;
	VectorXd fun_grad_z = gradOpZ * testfun;
	double tol = 1.e-10;
	bool passedx = fun_grad_x.isApproxToConstant(test_x, tol);
	bool passedy = fun_grad_y.isApproxToConstant(test_y, tol);
	bool passedz = fun_grad_z.isApproxToConstant(test_z, tol);
	std::cout << "Test of gradient x:" << (passedx ? " OK.": " Failed.") << std::endl;
	std::cout << "Test of gradient y:" << (passedy ? " OK.": " Failed.") << std::endl;
	std::cout << "Test of gradient z:" << (passedz ? " OK.": " Failed.") << std::endl;
	if (!passedx) {
		std::cout << "\nGradient x:\nmin : " << fun_grad_x.minCoeff() << std::endl;
		std::cout << "max : " << fun_grad_x.maxCoeff() << std::endl;
		std::cout << "Expected: " << test_x << "\n";
	}
	if (!passedy) {
		std::cout << "\nGradient y:\nmin : " << fun_grad_y.minCoeff() << std::endl;
		std::cout << "max : "	<< fun_grad_y.maxCoeff() << std::endl;
		std::cout << "Expected: " << test_y << "\n";
	}
	if (!passedz) {
		std::cout << "\nGradient z:\nmin : " << fun_grad_z.minCoeff() << std::endl;
		std::cout << "max : " << fun_grad_z.maxCoeff()	<< std::endl;
		std::cout << "Expected: " << test_z << "\n";
	}
}


void PreprocessMesh::prepareDirichletMatrix() {
	MeshBoundary BEmsh(msh.xyz, msh.fel);
	BEmsh.extractFeatures();
	msh.boundaryNodes = BEmsh.getBoundaryNodes();
	// Keep values at boundary nodes fixed.
	msh.dirichletMatrix = msh.stiff;
	VectorXi myDirichletNodes = Map<VectorXi>(msh.boundaryNodes.data(), msh.boundaryNodes.size());
	int bnx = msh.boundaryNodes.size();
	int nx = msh.dirichletMatrix.rows();
	if (msh.dirichletMatrix.IsRowMajor) {
/*
 * Row-major version:
*/
		SpMat_RM tL = msh.dirichletMatrix;
		for (int di = 0; di < bnx; ++di) {
				int fixedNode = myDirichletNodes(di);
				tL.row(fixedNode) = SpMat(nx, 1);
				tL.insert(fixedNode, fixedNode) = 1.;
//	        msh.dirichletMatrix.prune([fixedNode](int i, int, double) { return (i != fixedNode); }); // one-liner to set row(fixedNode) zero // This is terribly slow!
//		    msh.dirichletMatrix.insert(fixedNode, fixedNode) = 1.;
			}
			msh.dirichletMatrix = tL;
	} else {
/*
 * Column-major version:
 */
//	 	the matrix needs to be transposed. The .col() member has write access, the .row() member not.
		SpMat_CM tL = msh.dirichletMatrix.transpose();
		for (int di = 0; di < bnx; ++di) {
			int fixedNode = myDirichletNodes(di);
			tL.col(fixedNode) = SpMat_CM(nx, 1);
			tL.insert(fixedNode, fixedNode) = 1.;
//		msh.dirichletMatrix.prune([fixedNode](int i, int j, double) { return (i != fixedNode); }); // one-liner to set row(fixedNode) zero
//		msh.dirichletMatrix.insert(fixedNode, fixedNode) = 1.;
		}
		msh.dirichletMatrix = tL.transpose();
	}
	msh.dirichletMatrix.prune(0,0);
	msh.dirichletMatrix.makeCompressed();

}

