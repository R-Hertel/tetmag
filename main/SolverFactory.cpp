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
 * SolverFactory.cpp
 *
 *  Created on: Mar 7, 2017
 *      Author: riccardo
 */

#include "SolverFactory.h"
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/SparseLU>
#include <iostream>
#include <string>
#ifdef USE_CUDA
#include "GPUSolver.h"
#endif
#include "AMGCLSolver.h"
// #include "ViennaSolver.h"

using namespace Eigen;

void SolverFactory::setMatrix(const SpMat& A_ ) {
	A = A_; // columnMajor version still needed for GPU solver
	Ar = A_;
}


void SolverFactory::setLoadVector(const Ref<const VectorXd>& b_) {
	b = b_;
}


void SolverFactory::setCGTolerance(const double cgTol_) {
	cgTol = cgTol_;
}


void SolverFactory::definePreconditioner(std::string preconditioner) {
	preconditionerType = preconditioner;
}


VectorXd SolverFactory::result() {
	return x;
}


struct SparseStorage {
	VectorXi Index;
	VectorXi Ptr;
	VectorXd Values;
};


class LUSolver : public SolverFactory {
public:
//	SparseLU<SpMat_CM, COLAMDOrdering<int> > solver;
	SparseLU<SpMat_CM> solver;
	bool compute() {
		Ac = A;
		solver.compute(Ac);
		return true;
	}
	void solve() { x = solver.solve(b); }
	bool wasSuccessful() {
		return (solver.info() == Success); }
	void report() {
		std::cout << "solver " << ((solver.info() == Success) ? "successful." : "FAILED.") << std::endl;
	}
};


class CGSolver : public SolverFactory {
public:
//	BiCGSTAB<SpMat_RM, DiagonalPreconditioner<double> > solver;
	BiCGSTAB<SpMat_RM, IncompleteLUT<double> > solver;
//	BiCGSTAB<SpMat_RM > solver;
	bool compute() {
		Ar = A; // setting row-major matrix
		previous = VectorXd::Ones(Ar.rows()); // initial guess
		solver.setTolerance( cgTol );
//		std::cout << "CG solver tolerance: " << std::scientific << solver.tolerance() << std::fixed << std::endl;
		std::cout << "Preconditioning..." << std::flush;
		solver.compute(Ar);
		std::cout << " done." << std::endl;
		return true; // only to match factory signature
	}
	void solve(){
		x = solver.solveWithGuess(b, previous);
		previous = x;
	}
	void report() {
			std::cout << "#of iterations: " << solver.iterations() << std::endl;
			std::cout << "estimated error: " << solver.error()     << std::endl;
			std::cout << "solver " << ((solver.info() == Success) ? "successful." : "FAILED.") << std::endl;
	}
	bool wasSuccessful() { return (solver.info() == Success); }
};


std::shared_ptr<SolverFactory> SolverFactory::makeSolver(std::string type) {
	if (type == "lu") {
		return std::shared_ptr<SolverFactory>(std::make_shared<LUSolver>());
	} else if (( type == "cg" ) || ( type == "cpu" )) {
		return std::shared_ptr<SolverFactory>(std::make_shared<AMGCLSolver>());
	} else if ( type == "bcg"  ) {
		return std::shared_ptr<SolverFactory>(std::make_shared<CGSolver>());
	} else if (type == "gpu") {
#ifdef USE_CUDA
// this is the GPU version of the AMGCL solver:
		return std::shared_ptr<SolverFactory>(std::make_shared<GPUSolver>());
#else
		std::cout << "GPU option is not available in this version." << std::endl;
		return std::shared_ptr<SolverFactory>(std::make_shared<AMGCLSolver>());
#endif
	} else {
		std::cout << "solver option not understood." << std::endl;
		exit(0);
	}
}
