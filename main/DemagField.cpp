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
#include <vector>
#include <Eigen/SparseCore>
#include <Eigen/Dense>
#include <iostream>
#include "MeshData.h"
#ifdef USE_CUDA
#include "GpuPotential.h"
#endif
#include "h2interface.h"
#include "PhysicalConstants.h"

/*
#include <chrono>
using milli = std::chrono::milliseconds;
*/

using namespace Eigen;
enum coords {x, y, z};


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
		boundaryNodes(msh.boundaryNodes), dirichletBEM(msh.laplaceBEM), dirichletMatrix(msh.dirichletMatrix),
		neumannMatrix(msh.stiff), // only for initialization, actual matrix is set up in prepareNeumannMatrix()
		tGradX(msh.tGradX), tGradY(msh.tGradY), tGradZ(msh.tGradZ),
		gradX(msh.gradX), gradY(msh.gradY), gradZ(msh.gradZ),
		Js(Js_), NodeVolume(msh.NodeVolume),
		useH2(useH2_) {
	u1Timer.reset(); h2Timer.reset(); u2Timer.reset(), demagTimer.reset();
	nx = neumannMatrix.rows();
	Hdem = MatrixXd::Zero(nx,3);
	if (!useH2) {
		bnx = dirichletBEM.rows();
	} else {
		bnx = getNumberOfVertices();
	}
	fixedNeumannNode = static_cast<int>(nx / 2); // select arbitrary node which will remain fixed
	prepareNeumannMatrix();
	dirichletMatrix.makeCompressed();

	if (useH2) {
		selectedMVPtype = [this](VectorXd& v) -> VectorXd { return h2MVP(v); };
	} else {
		selectedMVPtype = [this](VectorXd& v) -> VectorXd { return denseMVP(v); };
	}
}


MatrixXd DemagField::calcField(MRef & mag) {
	return calcGradientField(calcPotential(mag));
}


double DemagField::getDemagEnergy(MRef & mag) {
	return ( -Js.cwiseProduct(NodeVolume).transpose() * Hdem.cwiseProduct(mag) ).sum() / 2.;
}


MatrixXd DemagField::calcGradientField(const Ref<const VectorXd>& u) {
	if (useGPU) {
#ifdef USE_CUDA
		Hdem = gpuField.calcGradientField(u);
#endif
	} else {
		Hdem.col(x) = -gradX * u ;
		Hdem.col(y) = -gradY * u ;
		Hdem.col(z) = -gradZ * u ;
	}
	return Hdem ;
}


VectorXd DemagField::getDivM(MRef & mag) { // volume charges, for graphics output
	VectorXd divM = ( gradX * mag.col(x)
					+ gradY * mag.col(y)
					+ gradZ * mag.col(z));
	return divM;
}


VectorXd DemagField::rhsPoisson(MRef & mag) { // here mag is \mu0*Ms
	VectorXd divM;
	if (useGPU) {
#ifdef USE_CUDA
		divM = gpuField.calcDivM(mag);
#endif
	} else {
		divM = (  tGradX * mag.col(x).cwiseProduct(Js)
				+ tGradY * mag.col(y).cwiseProduct(Js)
				+ tGradZ * mag.col(z).cwiseProduct(Js));
	}
	return divM / PhysicalConstants::mu0;
}


VectorXd DemagField::h2MVP(VectorXd& v) {
	pvector res = H2_mvp_sub(v.data(), bnx);
	return Map<VectorXd>(res, bnx);
}


VectorXd DemagField::denseMVP(VectorXd& v) {
	return dirichletBEM * v;
}


VectorXd DemagField::mvp(VectorXd& v) {
	return selectedMVPtype(v);
}


void DemagField::outputTimer() {
        std::cout << "u1 time [s]:\t" << u1Timer.durationInMus() / 1.e6 << "\t(" << u1Timer.durationInMus() / demagTimer.durationInMus() * 100. <<" %)"<< std::endl;
        std::cout << "h2 time [s]:\t" << h2Timer.durationInMus() / 1.e6 << "\t(" << h2Timer.durationInMus() / demagTimer.durationInMus() * 100. <<" %)" << std::endl;
        std::cout << "u2 time [s]:\t" << u2Timer.durationInMus() / 1.e6 << "\t(" << u2Timer.durationInMus() / demagTimer.durationInMus() * 100. <<" %)"<< std::endl;
        std::cout << "total time [s]:\t" <<  demagTimer.durationInMus() / 1.e6 << std::endl;
}


VectorXd DemagField::calcPotential(MRef & mag ) {
	 demagTimer.start();
	VectorXd u1, u2;

	    u1Timer.start();
		u1 = solvePoisson( rhsPoisson(mag) );
		u1Timer.add();

		h2Timer.start();
		VectorXd boundaryValues = boundaryIntegral(u1);
		h2Timer.add();

		u2Timer.start();
		u2 = solveLaplace( boundaryValues );
		u2Timer.add();
		demagTimer.add();
//	}
	return (u1 + u2);
}


VectorXd DemagField::boundaryIntegral(const Eigen::Ref<const VectorXd>& u1) {
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


VectorXd DemagField::solveLaplace(const Ref<const VectorXd>& boundaryValues ) {
	dirichletSolver->setLoadVector(boundaryValues);
	dirichletSolver->solve();
//	assert(dirichletSolver->wasSuccessful());
//	dirichletSolver->report();
//	VectorXd u2 = dirichletSolver->result();
	return dirichletSolver->result();
}


VectorXd DemagField::solvePoisson(const Ref<const VectorXd>& divM ) {
	neumannSolver->setLoadVector(divM);
	neumannSolver->solve();
//	assert (neumannSolver->wasSuccessful());
//	neumannSolver->report();
//	VectorXd u1 = neumannSolver->result();
	return neumannSolver->result();
}


void DemagField::initializeSolvers(std::string solverType, std::string preconditionerType, double cgTol_) {
	cgTol = cgTol_;
	if (solverType == "gpu" || solverType == "pl") {
		useGPU = true;
#ifdef USE_CUDA
		gpuField.setGradientMatrices(tGradX, tGradY, tGradZ, gradX, gradY, gradZ, Js);
#endif
	} else {
		useGPU = false;
	}

		dirichletSolver = SolverFactory::makeSolver( solverType );
		dirichletSolver->definePreconditioner( preconditionerType );
		dirichletSolver->setCGTolerance( cgTol );
		assert(dirichletMatrix.isCompressed());
		dirichletSolver->setMatrix(dirichletMatrix);
		bool success;
		success = dirichletSolver->compute();
		if (!success) {
			std::cout << "setup of Laplace solver failed." << std::endl;
		}

		neumannSolver = SolverFactory::makeSolver( solverType );
		neumannSolver->definePreconditioner( preconditionerType );
		neumannSolver->setCGTolerance( cgTol );
		assert(neumannMatrix.isCompressed());
		neumannSolver->setMatrix(neumannMatrix);
		success = neumannSolver->compute();
		if (!success) {
			std::cout << "setup of Poisson solver failed." << std::endl;
		}

		if (!success) {
			std::cout << "The selected solver option doesn't work for this problem.\n"
				"Please try setting the option\n\tsolver type = LU \n or\n\tsolver type = BCG" << std::endl;
			exit(0);
		}

//	}
}
