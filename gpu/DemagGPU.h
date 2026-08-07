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
 * DemagGPU.h
 *
 *  Created on: Aug 7, 2026
 *      Author: riccardo
 */

#ifndef GPU_DEMAGGPU_H_
#define GPU_DEMAGGPU_H_
#include <thrust/device_vector.h>
#include <Eigen/Dense>
#include <Eigen/SparseCore>
#include <memory>
#include "DemagBackend.h"
#include "SolverFactory.h"
#include "SpMatCUDA.h"
#include "typedefs.h"

class DemagGPU {
private:
	size_t nx;

	// -grad, so that Hdem follows directly from the potential:
	std::shared_ptr<SpMatCUDA> negGradX_cuda, negGradY_cuda, negGradZ_cuda;
	std::shared_ptr<SpMatCUDA> tGradX_cuda, tGradY_cuda, tGradZ_cuda;

	std::shared_ptr<devVecD> Js_dev;
	std::shared_ptr<devVecD> mx_dev, my_dev, mz_dev;
	std::shared_ptr<devVecD> Jx, Jy, Jz; // magnetic polarization Js*m at the nodes
	std::shared_ptr<devVecD> x_tmp, y_tmp, z_tmp;
	std::shared_ptr<devVecD> divM_d, u_dev;

	Eigen::VectorXd divM;
	Eigen::VectorXd u1, u2;
	Eigen::VectorXd boundaryValues;
	Eigen::MatrixXd Hdem;
	std::shared_ptr<SolverFactory> dirichletSolver;
	std::shared_ptr<SolverFactory> neumannSolver;
public:
	DemagGPU(size_t, const SpMat&, const SpMat&, const SpMat&,
			const SpMat&, const SpMat&, const SpMat&, const Eigen::VectorXd&);
	void setSolvers(std::shared_ptr<SolverFactory>, std::shared_ptr<SolverFactory>);
	void setMagnetization(MRef&);
	void computeRhs();
	void solvePoisson();
	const Eigen::VectorXd& poissonPotential();
	void setBoundaryValues(const Eigen::Ref<const Eigen::VectorXd>&);
	void solveLaplace();
	void computeField();
	const Eigen::MatrixXd& field();
};

std::shared_ptr<DemagBackend> makeDemagGPU(std::shared_ptr<DemagGPU>);

#endif /* GPU_DEMAGGPU_H_ */
