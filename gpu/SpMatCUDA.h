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
 * SpMatCUDA.h
 *
 *  Created on: Sep 17, 2020
 *      Author: hertel
 */

#ifndef GPU_SPMATCUDA_H_
#define GPU_SPMATCUDA_H_

#include <cusparse.h>
#include <Eigen/Dense>
#include <Eigen/SparseCore>
#include <thrust/device_vector.h>
#include <thrust/host_vector.h>
#include <vector>
#include "typedefs.h"

typedef thrust::host_vector<double> hostVecD;
typedef thrust::device_vector<double> devVecD;
typedef thrust::device_vector<int> devVecI;

devVecD devVecXd(const Eigen::VectorXd& );

class SpMatCUDA {
private:
	devVecD csrVals_d;
	devVecI csrRowPtr_d;
	devVecI csrColInd_d;
	cusparseDnVecDescr_t vecX, vecY;
	cusparseHandle_t handle;
	cusparseSpMatDescr_t matA;
	void *dBuffer;
	size_t bufferSize;
	int nnz, cols, rows;
	void setOnDev();
	void checkStatusCusparse(cusparseStatus_t&);
	template<class deviceType>
		void delete_vec(deviceType&);
public:
	const double alpha; // SAXPY parameters
	const double beta;
	SpMatCUDA(const SpMat&);
	void mvp(const devVecD&, devVecD&);
	SpMatCUDA();
	virtual ~SpMatCUDA();
};

#endif /* GPU_SPMATCUDA_H_ */
