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
 * DemagField.cpp
 *
 *  Created on: May 4, 2017
 *      Author: riccardo
 */

#include "DemagField.h"
#include "SolverFactory.h"
#include "DemagCPU.h"
#ifdef USE_CUDA
#include "DemagGPU.h"
#endif
#include <vector>
#include <cassert>
#include <Eigen/SparseCore>
#include <Eigen/Dense>
#include <iostream>
#include "MeshData.h"

/*
#include <chrono>
using milli = std::chrono::milliseconds;
*/

using namespace Eigen;


////////////////////////////////////////////////////////////////////////////////////////////////////
void checkIfSingular(SpMat stiff) {
	MatrixXd Mdense = MatrixXd( stiff );  // convert into dense matrix to calculate determinant
	std::cout << "Calculating determinant. Please wait..." << std::endl;
	double Mdet = Mdense.determinant();
	double small = 1.e-14;
	if (std::abs(Mdet) < small) {std::cout << "The matrix is SINGULAR." << std::endl;}
	std::cout << "Determinant: " << Mdet << std::endl;
}


DemagField::DemagField(const MeshData& msh, const VectorXd& Js_, bool useH2_ ) :
		dirichletMatrix(msh.dirichletMatrix),
		neumannMatrix(msh.stiff), // only for initialization, actual matrix is set up in prepareNeumannMatrix()
		tGradX(msh.tGradX), tGradY(msh.tGradY), tGradZ(msh.tGradZ),
		gradX(msh.gradX), gradY(msh.gradY), gradZ(msh.gradZ),
		Js(Js_), NodeVolume(msh.NodeVolume) {
	u1Timer.reset(); h2Timer.reset(); u2Timer.reset(), demagTimer.reset();
	nx = neumannMatrix.rows();
	bem = std::make_shared<BEMOperator>(msh, nx, useH2_);
	fixedNeumannNode = static_cast<int>(nx / 2); // select arbitrary node which will remain fixed
	prepareNeumannMatrix();
	dirichletMatrix.makeCompressed();
}


MatrixXd DemagField::calcField(MRef & mag) {
	demagTimer.start();

	u1Timer.start();
	backend->setMagnetization(mag);
	backend->computeRhs();
	backend->solvePoisson();
	u1Timer.add();

	h2Timer.start();
	backend->setBoundaryValues( bem->boundaryIntegral( backend->poissonPotential() ) );
	h2Timer.add();

	u2Timer.start();
	backend->solveLaplace();
	u2Timer.add();
	demagTimer.add();

	backend->computeField();
	return backend->field();
}


double DemagField::getDemagEnergy(MRef & mag) {
	return ( -Js.cwiseProduct(NodeVolume).transpose() * backend->field().cwiseProduct(mag) ).sum() / 2.;
}


void DemagField::outputTimer() {
        std::cout << "u1 time [s]:\t" << u1Timer.durationInMus() / 1.e6 << "\t(" << u1Timer.durationInMus() / demagTimer.durationInMus() * 100. <<" %)"<< std::endl;
        std::cout << "h2 time [s]:\t" << h2Timer.durationInMus() / 1.e6 << "\t(" << h2Timer.durationInMus() / demagTimer.durationInMus() * 100. <<" %)" << std::endl;
        std::cout << "u2 time [s]:\t" << u2Timer.durationInMus() / 1.e6 << "\t(" << u2Timer.durationInMus() / demagTimer.durationInMus() * 100. <<" %)"<< std::endl;
        std::cout << "total time [s]:\t" <<  demagTimer.durationInMus() / 1.e6 << std::endl;
}


void DemagField::prepareNeumannMatrix() {
	int k = fixedNeumannNode;
// set row(k) to zero:
	if (neumannMatrix.IsRowMajor) {  // Row-Major version
		SpMat_RM tNr = neumannMatrix;
		tNr.row(k) = SpMat(nx, 1);
		neumannMatrix = tNr;
	} else { // Column-Major version
		SpMat_CM tNc = neumannMatrix.transpose(); // required because of column-major storage and absence of .row() access
		tNc.col(k) = SpMat(nx, 1);
		neumannMatrix = tNc.transpose();
	}

	neumannMatrix.insert(k, k) = 1.;
	neumannMatrix.prune(0, 0);
	neumannMatrix.makeCompressed();
}


std::shared_ptr<SolverFactory> DemagField::makeAndSetupSolver(const std::string& solverType,
		const std::string& preconditionerType, const SpMat& matrix, const std::string& name) {
	std::shared_ptr<SolverFactory> solver = SolverFactory::makeSolver( solverType );
	solver->definePreconditioner( preconditionerType );
	solver->setCGTolerance( cgTol );
	assert(matrix.isCompressed());
	solver->setMatrix(matrix);
	if (!solver->compute()) {
		std::cerr << "setup of " << name << " solver failed.\n"
			"The selected solver option doesn't work for this problem.\n"
			"Please try setting the option\n\tsolver type = LU \n or\n\tsolver type = BCG" << std::endl;
		exit(1);
	}
	return solver;
}


void DemagField::initializeSolvers(std::string solverType, std::string preconditionerType, double cgTol_) {
	cgTol = cgTol_;
	std::shared_ptr<SolverFactory> dirichletSolver =
			makeAndSetupSolver(solverType, preconditionerType, dirichletMatrix, "Laplace");
	std::shared_ptr<SolverFactory> neumannSolver =
			makeAndSetupSolver(solverType, preconditionerType, neumannMatrix, "Poisson");

#ifdef USE_CUDA
	if (solverType == "gpu") {
		std::shared_ptr<DemagGPU> gpuImpl =
				std::make_shared<DemagGPU>(nx, tGradX, tGradY, tGradZ, gradX, gradY, gradZ, Js);
		gpuImpl->setSolvers(neumannSolver, dirichletSolver);
		backend = makeDemagGPU(gpuImpl);
		return;
	}
#endif
	std::shared_ptr<DemagCPU> cpuImpl =
			std::make_shared<DemagCPU>(nx, tGradX, tGradY, tGradZ, gradX, gradY, gradZ, Js);
	cpuImpl->setSolvers(neumannSolver, dirichletSolver);
	backend = makeDemagCPU(cpuImpl);
}
