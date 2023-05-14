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
 * GpuPotential.h
 *
 *  Created on: Jun 15, 2018
 *      Author: riccardo
 */

#ifndef GPU_GPUPOTENTIAL_H_
#define GPU_GPUPOTENTIAL_H_

#include <thrust/device_vector.h>
#include <thrust/host_vector.h>
#include <Eigen/Dense>
#include <Eigen/SparseCore>
#include "typedefs.h"
#include "SpMatCUDA.h"
#include <memory>

typedef double value_type;
typedef std::vector<value_type> host_vec;
typedef thrust::host_vector< value_type > vec_h;
typedef thrust::device_vector< value_type > dev_vec;
typedef thrust::device_vector<double> devVecD;
typedef thrust::device_vector<int> devVecI;

class GpuPotential{
private:

	std::vector<double> divM_h;
	std::vector<double> Hdem_vec;
	int nx;

	std::shared_ptr<SpMatCUDA> gradX_cuda;
	std::shared_ptr<SpMatCUDA> gradY_cuda;
	std::shared_ptr<SpMatCUDA> gradZ_cuda;

	std::shared_ptr<SpMatCUDA> tGradX_cuda;
	std::shared_ptr<SpMatCUDA> tGradY_cuda;
	std::shared_ptr<SpMatCUDA> tGradZ_cuda;

	std::shared_ptr<devVecD> Jx, Jy, Jz;
	std::shared_ptr<devVecD> divM_d, u_dev;
	std::shared_ptr<devVecD> Js_dev;
	std::shared_ptr<devVecD>  mx_dev, my_dev, mz_dev;
	std::shared_ptr<devVecD>  x_tmp, y_tmp, z_tmp;

public:
	void setGradientMatrices(const SpMat&, const SpMat&, const SpMat&, const SpMat&, const SpMat&, const SpMat&, const Eigen::VectorXd&);
	Eigen::VectorXd calcDivM(const Eigen::MatrixXd& );
	Eigen::Matrix<value_type, Eigen::Dynamic, Eigen::Dynamic> calcGradientField(const Eigen::VectorXd&);
};

#endif /* GPU_GPUPOTENTIAL_H_ */
