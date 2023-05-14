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
 * AMGCLSolver.cpp
 *
 *  Created on: Sep 4, 2020
 *      Author: hertel
 */

#include "AMGCLSolver.h"
#include <Eigen/Dense>
#include <cmath>
#include <Eigen/../unsupported/Eigen/SparseExtra>
using namespace Eigen;


void AMGCLSolver::init() {
      prm.precond.coarsening.relax = 0.75;
      prm.solver.tol = SolverFactory::cgTol;
      prm.solver.maxiter = 100;
// not sure about these parameters:
      prm.precond.coarse_enough = 100;
      prm.precond.direct_coarse = false;
//      amgcl::amg::params genpam;
//      genpam.direct_coarse = true;
  // for schur version:  
  // prm.precond.usolver.solver.tol = 1e-2;
  // prm.precond.psolver.solver.tol = 1e-2;
  // prm.precond.approx_schur = false;
//    std::cout << "CG solver tolerance: " << std::scientific << prm.solver.tol  << std::fixed << std::endl;
}

void AMGCLSolver::report() {
  std::cout << *solver << std::endl;
}

bool AMGCLSolver::setup(SpMat& A_) {
	bool success = true;
    SpMat_RM A = A_; // copy into / enforce RowMajor format
    int nnz = A.nonZeros();
    n = A.cols();
    b_n.resize(n);
    x_n.resize(n);
    for (int i = 0; i < n; ++i){
    	x_n[i] = 0.;
    }
    A.makeCompressed();
    std::vector<int> col(nnz);
    std::vector<int> ptr(A.outerSize() + 1);
    std::vector<double> val(nnz);

    x.resize(n);
//    x_h.resize(n);
//    b_h.resize(n);
    std::copy( A.valuePtr(), A.valuePtr( ) + nnz , val.begin());
    std::copy( A.innerIndexPtr(), A.innerIndexPtr() + nnz, col.begin() );
    std::copy( A.outerIndexPtr(), A.outerIndexPtr() + A.outerSize() + 1 , ptr.begin() );

    try{
    	solver = std::make_shared<AMGSolver>(std::tie(n,  ptr, col,  val), prm, bprm);
    } catch (...){
    	success = false;
    }
    if (!success) {
    	saveMarket(A, "solverMatrix.mtx");
    	std::cout << "Failed to set up matrix for solver." << std::endl;
    	exit(1);
    }
//    saveMarket(A, "LaplaceMatrix.mtx");
//    exit(0);
    return success;
}

// Skeletons
void AMGCLSolver::solve() {
	std::copy(b.data(), b.data() + n, b_n.data());
	std::copy(x.data(), x.data() + n, x_n.data());
	std::tie(iters, error) = solver->operator()(b_n, x_n);
//	std::cout << "iters = " << iters << ", err = " << error << "\n\n";
	std::copy(x_n.data(), x_n.data() + n, x.data());
}

bool AMGCLSolver::compute() {
	init();
	return setup(A);
};


bool AMGCLSolver::wasSuccessful(){
	return std::isfinite(error);
};

