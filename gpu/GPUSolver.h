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
 * GPUSolver.h
 *
 *  Created on: Sep 2, 2020
 *      Author: hertel
 */

#ifndef GPU_GPUSOLVER_H_
#define GPU_GPUSOLVER_H_
#include <amgcl/backend/cuda.hpp>
#include <amgcl/adapter/crs_tuple.hpp>
#include <amgcl/make_solver.hpp>
#include <amgcl/amg.hpp>
#include <amgcl/coarsening/smoothed_aggregation.hpp>
//#include <amgcl/relaxation/spai0.hpp>
#include <amgcl/relaxation/ilu0.hpp>
//#include <amgcl/profiler.hpp>
//#include <amgcl/solver/bicgstab.hpp>
#include <amgcl/solver/gmres.hpp>

#include <Eigen/SparseCore>
#include <vector>
#include <tuple>
#include <thrust/device_vector.h>
#include "SolverFactory.h"
#include "typedefs.h"


typedef amgcl::backend::cuda<double> Backend;

typedef amgcl::make_solver<
           amgcl::amg< Backend,
           amgcl::coarsening::smoothed_aggregation,
//         amgcl::relaxation::spai0
           amgcl::relaxation::ilu0
          >,
			amgcl::solver::gmres<Backend>
//          amgcl::solver::bicgstab<Backend>
        > Solver;

class GPUSolver : public SolverFactory {

public:
	GPUSolver();
	virtual ~GPUSolver();
private:
	  int n;
	  Backend::params bprm;
	  Solver::params prm;
	  thrust::device_vector<double> x_d;
	  thrust::device_vector<double> b_d;

	  int iters;
	  double error;
public:
	  void solve();
	  bool compute();
	  void report();
	  bool wasSuccessful();
	  bool setup(SpMat&);
	  void init();
	  std::shared_ptr<Solver> solver;
};


#endif /* GPU_GPUSOLVER_H_ */
