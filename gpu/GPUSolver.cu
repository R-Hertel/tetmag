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
 * GPUSolver.cu
 *
 *  Created on: Sep 2, 2020
 *      Author: hertel
 */

#include "GPUSolver.h"
#include <Eigen/Dense>
#include <Eigen/SparseCore>
#include <iostream>
#include <cmath>
GPUSolver::GPUSolver() {}


GPUSolver::~GPUSolver() {}


void GPUSolver::init() {
    cusparseCreate(&bprm.cusparse_handle);
    prm.precond.coarsening.relax = 0.75;
    prm.solver.tol  = SolverFactory::cgTol;
    prm.solver.maxiter = 100;
//    std::cout << "CG solver tolerance: " << std::scientific << prm.solver.tol  << std::fixed << std::endl;
}


void GPUSolver::report() {
  std::cout << *solver << std::endl;
}


bool GPUSolver::setup(SpMat& A_) {
	bool success = true;
    SpMat_RM A = A_; // copy into / enforce RowMajor format
    A.makeCompressed();
    int nnz = A.nonZeros();
    std::vector<int> col(nnz);
    std::vector<int> ptr(A.outerSize() + 1);
    std::vector<double> val(nnz);
    n = A.cols();
    x.resize(n);
    x_d.resize(n);
    b_d.resize(n);

    std::copy( A.valuePtr(), A.valuePtr( ) + nnz , val.begin());
    std::copy( A.innerIndexPtr(), A.innerIndexPtr() + nnz, col.begin() );
    std::copy( A.outerIndexPtr(), A.outerIndexPtr() + A.outerSize() + 1 , ptr.begin() );
    try {
    solver = std::make_shared<Solver>(std::tie(n,  ptr, col,  val), prm, bprm);
    } catch (...) {
    	success = false;
    }
    return success;
}


void GPUSolver::solve() {
	thrust::copy(b.data(), b.data() + n, b_d.begin());
	std::tie(iters, error) = solver->operator()(b_d, x_d);
	thrust::copy(x_d.begin(), x_d.end(), x.data());
}


bool GPUSolver::compute() {
	init();
	return setup(A);
}


bool GPUSolver::wasSuccessful(){
	return std::isfinite(error);
}

