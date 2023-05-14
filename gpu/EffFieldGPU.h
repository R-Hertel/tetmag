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
 * EffFieldGPU.h
 *
 *  Created on: Sep 4, 2019
 *      Author: riccardo
 */

#ifndef GPU_EFFFIELDGPU_H_
#define GPU_EFFFIELDGPU_H_

#include <cusparse.h>
#include <thrust/device_vector.h>
#include <thrust/host_vector.h>
#include <Eigen/Dense>
#include <Eigen/SparseCore>
#include <vector>
#include "typedefs.h"
#include "SpMatCUDA.h"
#include "memory.h"
#include "Timer.h"

typedef double value_type;
typedef std::vector<value_type> host_vec;
typedef thrust::host_vector< value_type > vec_h;
typedef thrust::device_vector< value_type > dev_vec;
typedef thrust::device_vector<double> devVecD;
typedef thrust::device_vector<int> devVecI;


class EffFieldGPU {
private:
	int nx;
	Eigen::VectorXd Hxc_unrolled;
	std::vector<value_type> retVecLLG;
	Eigen::Matrix<value_type, Eigen::Dynamic, Eigen::Dynamic> retMatXd;
	std::shared_ptr<SpMatCUDA> xc_cuda;
	std::shared_ptr<SpMatCUDA> GradX_cuda, GradY_cuda, GradZ_cuda;
	std::shared_ptr<SpMatCUDA> tGradX_cuda, tGradY_cuda, tGradZ_cuda;
	std::shared_ptr<dev_vec> Hxcx_d;
	std::shared_ptr<dev_vec> Hxcy_d;
	std::shared_ptr<dev_vec> Hxcz_d;
	std::shared_ptr<dev_vec> dxdt;
	std::shared_ptr<dev_vec> mag3_n;
	std::shared_ptr<dev_vec> mx_d;
	std::shared_ptr<dev_vec> my_d;
	std::shared_ptr<dev_vec> mz_d;
	std::shared_ptr<dev_vec> heffx_d, heffy_d, heffz_d;
	std::shared_ptr<dev_vec> eta_jx_d, eta_jy_d, eta_jz_d;
	std::shared_ptr<dev_vec> curlM, cmx1, cmx2, cmy1, cmy2, cmz1, cmz2;
	std::shared_ptr<dev_vec> dmi_fac, dmi3;
	std::shared_ptr<dev_vec>  nv_surf_x, nv_surf_y, nv_surf_z;
	std::shared_ptr<dev_vec> surfTerm;
	std::shared_ptr<dev_vec> u_xx, u_xy, u_xz, u_yx, u_yy, u_yz, u_zx, u_zy, u_zz;
	std::shared_ptr<dev_vec> ju_xx, ju_xy, ju_xz, ju_yx, ju_yy, ju_yz, ju_zx, ju_zy, ju_zz;
	std::shared_ptr<dev_vec> u_term_stt, ustt_d;
	std::shared_ptr<dev_vec> kuAxis_x, kuAxis_y, kuAxis_z, kuVal;
	std::shared_ptr<dev_vec> hani_x, hani_y, hani_z;
	std::shared_ptr<dev_vec> kx_mx, ky_my, kz_mz;
	std::shared_ptr<dev_vec> tmp0, tmp1;
	Timer copyTimer;
	void calcCurlM();
	dev_vec cwiseProduct(const dev_vec& , const dev_vec & );


public:

//  [
	std::vector< value_type > ClassicLLG_hst(MRef&, const value_type);
	std::vector< value_type > LLG_noPrec_hst(MRef&, const value_type);
	std::vector< value_type > STT_term_LLG_hst(MRef&, const value_type, const value_type);
// ]

	void setMagDev(MRef&);
	void setMagDev(const std::vector<value_type>&);
	void setMagDev(const devVecD&);
	void setExchangeMatOnDev(const SpMat& );
	void setGradientMatsOnDev(const SpMat&, const SpMat&, const SpMat& );

	Eigen::Matrix< value_type, Eigen::Dynamic, Eigen::Dynamic > ExchangeFieldGPU();
	Eigen::Matrix< value_type, Eigen::Dynamic, Eigen::Dynamic > UniaxialAnisotropyField();
	Eigen::Matrix< value_type, Eigen::Dynamic, Eigen::Dynamic > CubicAnisotropyField();

	void setSTTDataOnDevice(const SpMat&, const SpMat&, const SpMat&,
				const Eigen::VectorXd&,
				const Eigen::VectorXd&,
				const Eigen::VectorXd&);	
	Eigen::Matrix<value_type, Eigen::Dynamic, Eigen::Dynamic> UTermSTT_GPU();

	std::shared_ptr<dev_vec> ClassicLLG_dev(MRef&, const value_type);
	std::shared_ptr<dev_vec> LLG_noPrec_dev(MRef&, const value_type);
	std::shared_ptr<dev_vec> STT_term_LLG_dev(MRef&, const value_type, const value_type);

	value_type MaxTorque(MRef&);
	void setUniaxialAnisotropy(const Eigen::MatrixXd&, const Eigen::VectorXd&);
	void setCubicAnisotropy(const std::vector<Eigen::Matrix3d>&, const Eigen::VectorXd&, const Eigen::VectorXd&);
	void setDMIdata(const Eigen::VectorXd&, const Eigen::VectorXd&, const Eigen::MatrixXd&, const Eigen::VectorXd&);
	Eigen::Matrix<value_type, Eigen::Dynamic, Eigen::Dynamic> DMIField();
	void NormalizeMag(Eigen::MatrixXd&, int);
	void init(int);

	EffFieldGPU(const int );
	void displayTimer();
};

#endif /* GPU_EFFFIELDGPU_H_ */
