/*
    tetmag - A general-purpose finite-element micromagnetic simulation software package
    Copyright (C) 2016-2026 CNRS and Université de Strasbourg

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
 * BEMOperator.cpp
 *
 *  Created on: Aug 7, 2026
 *      Author: riccardo
 */

#include "BEMOperator.h"
#include "MeshData.h"
#include "h2interface.h"
#include <Eigen/Dense>
#include <vector>

using namespace Eigen;


BEMOperator::BEMOperator(const MeshData& msh, size_t nx_, bool useH2) :
		nx(nx_), boundaryNodes(msh.boundaryNodes), dirichletBEM(msh.laplaceBEM) {
	if (!useH2) {
		bnx = dirichletBEM.rows();
		mvp = [this](VectorXd& v) -> VectorXd { return denseMVP(v); };
	} else {
		bnx = getNumberOfVertices();
		mvp = [this](VectorXd& v) -> VectorXd { return h2MVP(v); };
	}
}


VectorXd BEMOperator::h2MVP(VectorXd& v) {
	pvector res = H2_mvp_sub(v.data(), bnx);
	return Map<VectorXd>(res, bnx);
}


VectorXd BEMOperator::denseMVP(VectorXd& v) {
	return dirichletBEM * v;
}


VectorXd BEMOperator::boundaryIntegral(const Eigen::Ref<const VectorXd>& u1) {
	VectorXd boundaryValues = VectorXd::Zero(nx);
	VectorXd u1Boundary(bnx);
#ifdef _OPENMP
#pragma omp parallel for
#endif
	for (size_t i = 0; i < bnx; ++i) {
		u1Boundary(i) = u1(boundaryNodes[i]);
	}
	VectorXd u2Boundary = mvp(u1Boundary);
#ifdef _OPENMP
#pragma omp parallel for
#endif
	for (size_t i = 0; i < bnx; ++i) {
		boundaryValues(boundaryNodes[i]) = u2Boundary(i);
	}
	return boundaryValues;
}
