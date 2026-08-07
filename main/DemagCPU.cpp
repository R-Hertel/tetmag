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
 * DemagCPU.cpp
 *
 *  Created on: Aug 7, 2026
 *      Author: riccardo
 */

#include "DemagCPU.h"
#include "PhysicalConstants.h"
#include <Eigen/Dense>

using namespace Eigen;
enum coords {x, y, z};


DemagCPU::DemagCPU(size_t nx_, const SpMat& tGradX_, const SpMat& tGradY_, const SpMat& tGradZ_,
		const SpMat& gradX_, const SpMat& gradY_, const SpMat& gradZ_, const VectorXd& Js_) :
		nx(nx_), tGradX(tGradX_), tGradY(tGradY_), tGradZ(tGradZ_),
		negGradX(-gradX_), negGradY(-gradY_), negGradZ(-gradZ_), Js(Js_) {
	Jx.resize(nx);
	Jy.resize(nx);
	Jz.resize(nx);
	Hdem = MatrixXd::Zero(nx, 3);
}


void DemagCPU::setSolvers(std::shared_ptr<SolverFactory> neumann, std::shared_ptr<SolverFactory> dirichlet) {
	neumannSolver = neumann;
	dirichletSolver = dirichlet;
}


void DemagCPU::setMagnetization(MRef& mag) {
	Jx = mag.col(x).cwiseProduct(Js);
	Jy = mag.col(y).cwiseProduct(Js);
	Jz = mag.col(z).cwiseProduct(Js);
}


void DemagCPU::computeRhs() {
	divM = ( tGradX * Jx + tGradY * Jy + tGradZ * Jz ) / PhysicalConstants::mu0;
}


void DemagCPU::solvePoisson() {
	neumannSolver->setLoadVector(divM);
	neumannSolver->solve();
	u1 = neumannSolver->result();
}


const VectorXd& DemagCPU::poissonPotential() {
	return u1;
}


void DemagCPU::setBoundaryValues(const Eigen::Ref<const VectorXd>& boundaryValues_) {
	boundaryValues = boundaryValues_;
}


void DemagCPU::solveLaplace() {
	dirichletSolver->setLoadVector(boundaryValues);
	dirichletSolver->solve();
	u2 = dirichletSolver->result();
}


void DemagCPU::computeField() {
	const VectorXd u = u1 + u2;
	Hdem.col(x) = negGradX * u;
	Hdem.col(y) = negGradY * u;
	Hdem.col(z) = negGradZ * u;
}


const MatrixXd& DemagCPU::field() {
	return Hdem;
}


std::shared_ptr<DemagBackend> makeDemagCPU(std::shared_ptr<DemagCPU> impl) {
	std::shared_ptr<DemagBackend> ops = std::make_shared<DemagBackend>();
	ops->setMagnetization  = [impl](MRef& m) -> void { impl->setMagnetization(m); };
	ops->computeRhs        = [impl]() -> void { impl->computeRhs(); };
	ops->solvePoisson      = [impl]() -> void { impl->solvePoisson(); };
	ops->poissonPotential  = [impl]() -> const VectorXd& { return impl->poissonPotential(); };
	ops->setBoundaryValues = [impl](const Eigen::Ref<const VectorXd>& b) -> void { impl->setBoundaryValues(b); };
	ops->solveLaplace      = [impl]() -> void { impl->solveLaplace(); };
	ops->computeField      = [impl]() -> void { impl->computeField(); };
	ops->field             = [impl]() -> const MatrixXd& { return impl->field(); };
	return ops;
}
