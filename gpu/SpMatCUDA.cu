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
 * SpMatCUDA.cu
 *
 *  Created on: Sep 17, 2020
 *      Author: hertel
 */

#include "SpMatCUDA.h"
#include <iostream>
#include <vector>
#include <cuda_runtime_api.h>

using namespace Eigen;
// typedef SparseMatrix<double, ColMajor> SpMat_CM;
typedef thrust::device_vector<double> devVecD;

// construct thrust vector from Eigen Vector:
devVecD devVecXd(const Eigen::VectorXd& eig) {
   devVecD dv(std::vector<double>(eig.data(), eig.data() + eig.size()));
   return dv;
}


template<class deviceType>
void SpMatCUDA::delete_vec( deviceType& x_d ) {
	x_d.clear();
	x_d.shrink_to_fit();
	x_d.~device_vector();
}


SpMatCUDA::SpMatCUDA() : alpha(1.0), beta(0.0) {}


SpMatCUDA::SpMatCUDA( const SpMat& m1_) : alpha(1.0), beta(0.0) {
	SpMat_CM m1 = m1_; // copy into / enforce ColMajor format
	m1.makeCompressed();
	nnz = m1.nonZeros();
	cols = m1.cols();
	rows = m1.rows();
	Eigen::VectorXd cscValA_h    = Map<VectorXd>( m1.valuePtr(), nnz );
	Eigen::VectorXi	cscRowIndA_h = Map<VectorXi>( m1.innerIndexPtr(), nnz );
	Eigen::VectorXi cscColPtrA_h = Map<VectorXi>( m1.outerIndexPtr(), m1.outerSize() + 1 );
	cscVals_d.resize(cscValA_h.size());
	cscVals_d = devVecD ( cscValA_h.data(), cscValA_h.data() + cscValA_h.size() );
	cscCols_d.resize(cscColPtrA_h.size());
	cscCols_d = devVecI ( cscColPtrA_h.data(), cscColPtrA_h.data() + cscColPtrA_h.size() );
	cscRows_d.resize(cscRowIndA_h.size());
	cscRows_d = devVecI ( cscRowIndA_h.data(), cscRowIndA_h.data() + cscRowIndA_h.size() );
	setOnDev();
}


void SpMatCUDA::checkStatusCusparse(cusparseStatus_t& status) {
  if (status != CUSPARSE_STATUS_SUCCESS) {
    printf("CUSPARSE API failed at line %d with error: %s (%d)\n",
	   __LINE__, cusparseGetErrorString(status), status);
  } //else { std::cout << "CuSparseStatus = 0K" << std::endl;  }
}


void SpMatCUDA::setOnDev() {
  handle = NULL;
  cusparseCreate(&handle);
  cusparseStatus_t status;
  dBuffer = NULL;
  bufferSize = 0;
  VectorXd zero = VectorXd::Zero(rows);

  //	cusparseStatus_t status __attribute__ ((unused)) ; // this was in the old version.

  // prepare sparse matrix:
  status =  cusparseCreateCsr(&matA, rows, cols, nnz,
			      thrust::raw_pointer_cast(cscCols_d.data()),
			      thrust::raw_pointer_cast(cscRows_d.data()),
			      thrust::raw_pointer_cast(cscVals_d.data()),
			      CUSPARSE_INDEX_32I, CUSPARSE_INDEX_32I,
			      CUSPARSE_INDEX_BASE_ZERO, CUDA_R_64F);
  assert(status == CUSPARSE_STATUS_SUCCESS);
  checkStatusCusparse(status);

  // prepare output vector:
  cudaMalloc((void**) &dY, rows * sizeof(double));
  cudaMemcpy(dY, zero.data(), rows * sizeof(double), cudaMemcpyHostToDevice);
  cusparseCreateDnVec(&vecY, rows, dY, CUDA_R_64F);

  // prepare input vector:
  zero.resize(cols);
  cudaMalloc((void**) &dX, cols * sizeof(double));
  cudaMemcpy(dX, zero.data(), cols * sizeof(double), cudaMemcpyHostToDevice);
  cusparseCreateDnVec(&vecX, cols, dX, CUDA_R_64F);

  // prepare buffer:
  status = cusparseSpMV_bufferSize( handle, CUSPARSE_OPERATION_TRANSPOSE,
				   &alpha, matA, vecX, &beta, vecY, CUDA_R_64F,
#if (CUDART_VERSION > 11000)
				   CUSPARSE_SPMV_ALG_DEFAULT, &bufferSize) ;
#else
                                   CUSPARSE_MV_ALG_DEFAULT, &bufferSize) ;
#endif
  
  checkStatusCusparse(status);
  cudaMalloc(&dBuffer, bufferSize);
}


void SpMatCUDA::mvp() {
  cusparseStatus_t stat = cusparseSpMV(handle, CUSPARSE_OPERATION_TRANSPOSE,
				       &alpha, matA, vecX, &beta, vecY, CUDA_R_64F,
#if (CUDART_VERSION > 11000)
				   CUSPARSE_SPMV_ALG_DEFAULT, dBuffer) ;
#else
                                   CUSPARSE_MV_ALG_DEFAULT, dBuffer) ;
#endif
//  checkStatusCusparse(stat);
}


devVecD SpMatCUDA::mvp(const devVecD& x) {
	setX(x);
	mvp();
	return mvpResDev();
}


void SpMatCUDA::setX(const VectorXd& x) {
  cudaMemcpy(dX, x.data(), x.size() * sizeof(double),  cudaMemcpyHostToDevice);
}


void SpMatCUDA::setX(const devVecD& x) {
  cudaMemcpy(dX, thrust::raw_pointer_cast(x.data()), x.size() * sizeof(double),  cudaMemcpyDeviceToDevice);
}


VectorXd  SpMatCUDA::mvpResEig(){
  VectorXd output(rows);
   cudaMemcpy(output.data(), dY, rows * sizeof(double),  cudaMemcpyDeviceToHost);
   return output;
}


devVecD SpMatCUDA::mvpResDev(){
	return thrust::device_vector< double > (dY, dY + rows);
}


SpMatCUDA::~SpMatCUDA() {
	cusparseDestroySpMat (matA);
	cusparseDestroyDnVec (vecX);
	cusparseDestroyDnVec (vecY);
	cusparseDestroy (handle);
	delete_vec(cscVals_d);
	delete_vec(cscCols_d);
	delete_vec(cscRows_d);
}

