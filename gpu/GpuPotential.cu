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
 * GpuPotential.cu
 *
 *  Created on: Jun 15, 2018
 *      Author: riccardo
 */

#include "GpuPotential.h"
#include <iostream>
#include <Eigen/Dense>
#include "SpMatCUDA.h"
#include <thrust/transform.h>
#include <thrust/functional.h>

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

void GpuPotential::setGradientMatrices(const SpMat& tGradX, const SpMat& tGradY, const SpMat& tGradZ, 
									   const SpMat& gradX,  const SpMat& gradY,  const SpMat& gradZ, const VectorXd& Js) 
{
	gradX_cuda = std::make_shared<SpMatCUDA>( gradX );
	gradY_cuda = std::make_shared<SpMatCUDA>( gradY );
	gradZ_cuda = std::make_shared<SpMatCUDA>( gradZ );

	tGradX_cuda = std::make_shared<SpMatCUDA>( tGradX );
	tGradY_cuda = std::make_shared<SpMatCUDA>( tGradY );
	tGradZ_cuda = std::make_shared<SpMatCUDA>( tGradZ );

	nx = gradX.rows();

	mx_dev = std::make_shared<devVecD>(nx);
	my_dev = std::make_shared<devVecD>(nx);
	mz_dev = std::make_shared<devVecD>(nx);
	x_tmp  = std::make_shared<devVecD>(nx);
	y_tmp  = std::make_shared<devVecD>(nx);
	z_tmp  = std::make_shared<devVecD>(nx);
	divM_h.resize(nx);
	Hdem_vec.resize(3 * nx);

	divM_d = std::make_shared<devVecD>(nx);
	u_dev = std::make_shared<devVecD>(nx);
	Jx = std::make_shared<devVecD>(nx);
	Jy = std::make_shared<devVecD>(nx);
	Jz = std::make_shared<devVecD>(nx);
	Js_dev= std::make_shared<devVecD>(nx);
	thrust::copy(Js.data(), Js.data() + nx, Js_dev->begin());
}


Matrix<value_type, Dynamic, 1> GpuPotential::calcDivM(const MatrixXd& mag) {
	thrust::copy( mag.col(x).data(), mag.col(x).data() + nx, mx_dev->begin() );
	thrust::copy( mag.col(y).data(), mag.col(y).data() + nx, my_dev->begin() );
	thrust::copy( mag.col(z).data(), mag.col(z).data() + nx, mz_dev->begin() );

	thrust::transform(mx_dev->begin(), mx_dev->end(), Js_dev->begin(), Jx->begin(), thrust::multiplies<double>());
	thrust::transform(my_dev->begin(), my_dev->end(), Js_dev->begin(), Jy->begin(), thrust::multiplies<double>());
	thrust::transform(mz_dev->begin(), mz_dev->end(), Js_dev->begin(), Jz->begin(), thrust::multiplies<double>());

	*x_tmp = tGradX_cuda->mvp(*Jx);
	*y_tmp = tGradY_cuda->mvp(*Jy);
	*z_tmp = tGradZ_cuda->mvp(*Jz);

	thrust::transform(
			thrust::make_zip_iterator(thrust::make_tuple(x_tmp->begin(), y_tmp->begin(), z_tmp->begin() )),
			thrust::make_zip_iterator(thrust::make_tuple(x_tmp->end()  , y_tmp->end()  , z_tmp->end()   )),
			divM_d->begin(), sumOfThree<value_type>() );

	thrust::copy(divM_d->begin(), divM_d->end(), divM_h.begin());
	return Map<const Matrix<value_type, Dynamic, 1> >(divM_h.data(), nx); // If value_type is 'double', the return type is VectorXd
}


Matrix<value_type, Dynamic, Dynamic> GpuPotential::calcGradientField(const VectorXd& u){
	thrust::copy(u.data(), u.data() + nx, u_dev->begin());

	*x_tmp = gradX_cuda->mvp(*u_dev);
	*y_tmp = gradY_cuda->mvp(*u_dev);
	*z_tmp = gradZ_cuda->mvp(*u_dev);

	thrust::copy(x_tmp->begin(), x_tmp->end(), Hdem_vec.begin() );
	thrust::copy(y_tmp->begin(), y_tmp->end(), Hdem_vec.begin() + nx );
	thrust::copy(z_tmp->begin(), z_tmp->end(), Hdem_vec.begin() + 2 * nx );

	thrust::transform(Hdem_vec.begin(), Hdem_vec.end(), Hdem_vec.begin(), thrust::negate<value_type>());
	return Map<const Matrix<value_type, Dynamic, Dynamic> >(Hdem_vec.data(), nx, 3);
}

