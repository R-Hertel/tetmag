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
 * DemagGPU.cu
 *
 *  Created on: Aug 7, 2026
 *      Author: riccardo
 */

#include "DemagGPU.h"
#include "PhysicalConstants.h"
#include "SpMatCUDA.h"
#include <Eigen/Dense>
#include <thrust/copy.h>
#include <thrust/transform.h>
#include <thrust/functional.h>
#include <thrust/tuple.h>

using namespace Eigen;
enum coords {x, y, z};


template<class T>
struct sumOfThree {
    __host__ __device__
        T operator()(const thrust::tuple<T,T,T>& a) const
        {
            return thrust::get<0>(a)  +
                   thrust::get<1>(a)  +
                   thrust::get<2>(a) ;
        }
};


DemagGPU::DemagGPU(size_t nx_, const SpMat& tGradX, const SpMat& tGradY, const SpMat& tGradZ,
		const SpMat& gradX, const SpMat& gradY, const SpMat& gradZ, const VectorXd& Js) :
		nx(nx_) {
	negGradX_cuda = std::make_shared<SpMatCUDA>( -gradX );
	negGradY_cuda = std::make_shared<SpMatCUDA>( -gradY );
	negGradZ_cuda = std::make_shared<SpMatCUDA>( -gradZ );

	tGradX_cuda = std::make_shared<SpMatCUDA>( tGradX );
	tGradY_cuda = std::make_shared<SpMatCUDA>( tGradY );
	tGradZ_cuda = std::make_shared<SpMatCUDA>( tGradZ );

	mx_dev = std::make_shared<devVecD>(nx);
	my_dev = std::make_shared<devVecD>(nx);
	mz_dev = std::make_shared<devVecD>(nx);
	x_tmp  = std::make_shared<devVecD>(nx);
	y_tmp  = std::make_shared<devVecD>(nx);
	z_tmp  = std::make_shared<devVecD>(nx);

	divM_d = std::make_shared<devVecD>(nx);
	u_dev = std::make_shared<devVecD>(nx);
	Jx = std::make_shared<devVecD>(nx);
	Jy = std::make_shared<devVecD>(nx);
	Jz = std::make_shared<devVecD>(nx);
	Js_dev = std::make_shared<devVecD>(nx);
	thrust::copy(Js.data(), Js.data() + nx, Js_dev->begin());

	divM.resize(nx);
	Hdem = MatrixXd::Zero(nx, 3);
}


void DemagGPU::setSolvers(std::shared_ptr<SolverFactory> neumann, std::shared_ptr<SolverFactory> dirichlet) {
	neumannSolver = neumann;
	dirichletSolver = dirichlet;
}


void DemagGPU::setMagnetization(MRef& mag) {
	thrust::copy( mag.col(x).data(), mag.col(x).data() + nx, mx_dev->begin() );
	thrust::copy( mag.col(y).data(), mag.col(y).data() + nx, my_dev->begin() );
	thrust::copy( mag.col(z).data(), mag.col(z).data() + nx, mz_dev->begin() );

	thrust::transform(mx_dev->begin(), mx_dev->end(), Js_dev->begin(), Jx->begin(), thrust::multiplies<double>());
	thrust::transform(my_dev->begin(), my_dev->end(), Js_dev->begin(), Jy->begin(), thrust::multiplies<double>());
	thrust::transform(mz_dev->begin(), mz_dev->end(), Js_dev->begin(), Jz->begin(), thrust::multiplies<double>());
}


void DemagGPU::computeRhs() {
	tGradX_cuda->mvp(*Jx, *x_tmp);
	tGradY_cuda->mvp(*Jy, *y_tmp);
	tGradZ_cuda->mvp(*Jz, *z_tmp);

	thrust::transform(
			thrust::make_zip_iterator(thrust::make_tuple(x_tmp->begin(), y_tmp->begin(), z_tmp->begin() )),
			thrust::make_zip_iterator(thrust::make_tuple(x_tmp->end()  , y_tmp->end()  , z_tmp->end()   )),
			divM_d->begin(), sumOfThree<double>() );

	thrust::copy(divM_d->begin(), divM_d->end(), divM.data());
	divM /= PhysicalConstants::mu0;
}


void DemagGPU::solvePoisson() {
	neumannSolver->setLoadVector(divM);
	neumannSolver->solve();
	u1 = neumannSolver->result();
}


const VectorXd& DemagGPU::poissonPotential() {
	return u1;
}


void DemagGPU::setBoundaryValues(const Eigen::Ref<const VectorXd>& boundaryValues_) {
	boundaryValues = boundaryValues_;
}


void DemagGPU::solveLaplace() {
	dirichletSolver->setLoadVector(boundaryValues);
	dirichletSolver->solve();
	u2 = dirichletSolver->result();
}


void DemagGPU::computeField() {
	const VectorXd u = u1 + u2;
	thrust::copy(u.data(), u.data() + nx, u_dev->begin());

	negGradX_cuda->mvp(*u_dev, *x_tmp);
	negGradY_cuda->mvp(*u_dev, *y_tmp);
	negGradZ_cuda->mvp(*u_dev, *z_tmp);

	thrust::copy(x_tmp->begin(), x_tmp->end(), Hdem.col(x).data());
	thrust::copy(y_tmp->begin(), y_tmp->end(), Hdem.col(y).data());
	thrust::copy(z_tmp->begin(), z_tmp->end(), Hdem.col(z).data());
}


const MatrixXd& DemagGPU::field() {
	return Hdem;
}


std::shared_ptr<DemagBackend> makeDemagGPU(std::shared_ptr<DemagGPU> impl) {
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
