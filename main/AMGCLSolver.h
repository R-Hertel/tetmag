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
 * AMGCLSolver.h
 *
 *  Created on: Sep 4, 2020
 *      Author: hertel
 */

#ifndef MAIN_AMGCLSOLVER_H_
#define MAIN_AMGCLSOLVER_H_

#include "SolverFactory.h"

#include <amgcl/backend/builtin.hpp>
#include <amgcl/adapter/crs_tuple.hpp>
#include <amgcl/make_solver.hpp>
#include <amgcl/amg.hpp>
#include <amgcl/coarsening/smoothed_aggregation.hpp>
#include <amgcl/coarsening/smoothed_aggr_emin.hpp>
#include <amgcl/coarsening/aggregation.hpp>
#include <amgcl/preconditioner/schur_pressure_correction.hpp>
// #include <amgcl/relaxation/spai0.hpp>
#include <amgcl/relaxation/damped_jacobi.hpp>
#include <amgcl/relaxation/as_preconditioner.hpp>
#include <amgcl/relaxation/gauss_seidel.hpp>
#include <amgcl/relaxation/chebyshev.hpp>
#include <amgcl/relaxation/ilu0.hpp>
// #include <amgcl/relaxation/ilut.hpp>
#include <amgcl/profiler.hpp>
#include <amgcl/solver/cg.hpp>
#include <amgcl/solver/bicgstab.hpp>
#include <amgcl/solver/bicgstabl.hpp>
#include <amgcl/solver/gmres.hpp>
#include <amgcl/solver/fgmres.hpp>

#ifdef _OPENMP
#include <omp.h>
#endif

#include <Eigen/SparseCore>
#include <vector>
#include <tuple>

typedef amgcl::backend::builtin<double> CPU_Backend;
/*
// discussed here: https://github.com/ddemidov/amgcl/issues/144
typedef
amgcl::make_solver<
  amgcl::preconditioner::schur_pressure_correction<
  amgcl::make_solver<
    amgcl::amg<	CPU_Backend,
		amgcl::coarsening::smoothed_aggregation,
		amgcl::relaxation::ilu0
	>,
    amgcl::solver::cg<CPU_Backend>
	>,
  amgcl::make_solver<
    amgcl::relaxation::as_preconditioner<CPU_Backend, amgcl::relaxation::damped_jacobi>,
    amgcl::solver::cg<CPU_Backend>
	>
	>,
  amgcl::solver::fgmres<CPU_Backend>
  > AMGSolver;
*/  

typedef amgcl::make_solver<
  amgcl::amg<CPU_Backend,
	     //		amgcl::coarsening::aggregation,
	     amgcl::coarsening::smoothed_aggregation,
	     //		amgcl::coarsening::smoothed_aggr_emin,
	     //		amgcl::relaxation::spai0
	     //		amgcl::relaxation::gauss_seidel
	     //		amgcl::relaxation::damped_jacobi
	     //		amgcl::relaxation::chebyshev
	     amgcl::relaxation::ilu0
	     //	        amgcl::relaxation::ilut
	     >,
  amgcl::solver::gmres<CPU_Backend> > AMGSolver;
//		amgcl::solver::bicgstabl<CPU_Backend> > AMGSolver;
//		amgcl::solver::bicgstab<CPU_Backend> > AMGSolver;

  
class AMGCLSolver: public SolverFactory {
private:
		  int n;
		  CPU_Backend::params bprm;
		  AMGSolver::params prm;
//		  std::vector<double> x_h;
//	 	  std::vector<double> b_h;
		  amgcl::backend::numa_vector<double> b_n;
		  amgcl::backend::numa_vector<double> x_n;
		  int iters;
		  double error;
public:
		  void solve();
		  bool compute();
		  void report();
		  bool wasSuccessful();
		  bool setup(SpMat&);
		  void init();
		  std::shared_ptr<AMGSolver> solver;
};

#endif /* MAIN_AMGCLSOLVER_H_ */
