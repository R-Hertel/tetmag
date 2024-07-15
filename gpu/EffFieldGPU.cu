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
 * EffFieldGPU.cu
 *
 *  Created on: Sep 4, 2019
 *      Author: riccardo
 */

#include "typedefs.h"
#include "EffFieldGPU.h"
#include <cusparse.h>
#include <thrust/transform.h>
#include <thrust/transform_reduce.h>
#include <thrust/functional.h>
#include <thrust/fill.h>
#include <thrust/copy.h>
#include <thrust/iterator/zip_iterator.h>
#include <iostream>
using namespace Eigen;


enum Coords {
	x, y, z
};


EffFieldGPU::EffFieldGPU(const int nx_) : nx(nx_) {
	copyTimer.reset();
	init(nx);
}

void EffFieldGPU::setUniaxialAnisotropy(const MatrixXd& KuAxis_eig, const VectorXd& Ku1) {
	kuAxis_x = std::make_shared<dev_vec>(nx);
	kuAxis_y = std::make_shared<dev_vec>(nx);
	kuAxis_z = std::make_shared<dev_vec>(nx);
	kuVal = std::make_shared<dev_vec>(nx);
	thrust::copy(KuAxis_eig.col(x).data(), KuAxis_eig.col(x).data() + nx, kuAxis_x->begin());
	thrust::copy(KuAxis_eig.col(y).data(), KuAxis_eig.col(y).data() + nx, kuAxis_y->begin());
	thrust::copy(KuAxis_eig.col(z).data(), KuAxis_eig.col(z).data() + nx, kuAxis_z->begin());
	thrust::copy(Ku1.data(), Ku1.data() + nx , kuVal->begin());
}

void EffFieldGPU::setCubicAnisotropy(const std::vector<Matrix3d>& cubicAxes, const VectorXd& Kc1, const VectorXd& Kc2) {
// to be completed
}


Matrix<value_type, Dynamic, Dynamic> EffFieldGPU::ExchangeFieldGPU() {

	*Hxcx_d = xc_cuda->mvp(*mx_d);
	*Hxcy_d = xc_cuda->mvp(*my_d);
	*Hxcz_d = xc_cuda->mvp(*mz_d);

	copyTimer.start();
	thrust::copy( Hxcx_d->begin(), Hxcx_d->end() , Hxc_unrolled.data() );
  	thrust::copy( Hxcy_d->begin(), Hxcy_d->end() , Hxc_unrolled.data() + nx );
	thrust::copy( Hxcz_d->begin(), Hxcz_d->end() , Hxc_unrolled.data() + 2 * nx );
	copyTimer.end();
	copyTimer.add();

	return Map<Matrix<value_type, Dynamic, Dynamic> >( Hxc_unrolled.data(), nx, 3 );

}


template<class T>
struct sumOfThree {
    __host__ __device__
        T operator()(const thrust::tuple<T,T,T>& a) const
        {
            return thrust::get<0>(a)  +    // first vector
                   thrust::get<1>(a)  +    // second vector
                   thrust::get<2>(a) ;     // third vector
        }
};

dev_vec EffFieldGPU::cwiseProduct(const dev_vec& v1, const dev_vec& v2) {
  dev_vec out(v1.size());
  thrust::transform( v1.begin(), v1.end(), v2.begin(), out.begin(), thrust::multiplies<value_type>() );
  return out;
}


Matrix<double, Dynamic, Dynamic> EffFieldGPU::UTermSTT_GPU() {

  *u_xx = GradX_cuda->mvp(*mx_d);
  *u_xy = GradY_cuda->mvp(*mx_d);
  *u_xz = GradZ_cuda->mvp(*mx_d);

  *u_yx = GradX_cuda->mvp(*my_d);
  *u_yy = GradY_cuda->mvp(*my_d);
  *u_yz = GradZ_cuda->mvp(*my_d);

  *u_zx = GradX_cuda->mvp(*mz_d);
  *u_zy = GradY_cuda->mvp(*mz_d);
  *u_zz = GradZ_cuda->mvp(*mz_d);
  
  *ju_xx = cwiseProduct( *eta_jx_d, *u_xx );
  *ju_xy = cwiseProduct( *eta_jy_d, *u_xy );
  *ju_xz = cwiseProduct( *eta_jz_d, *u_xz );

  *ju_yx = cwiseProduct( *eta_jx_d, *u_yx );
  *ju_yy = cwiseProduct( *eta_jy_d, *u_yy );
  *ju_yz = cwiseProduct( *eta_jz_d, *u_yz );
  
  *ju_zx = cwiseProduct( *eta_jx_d, *u_zx );
  *ju_zy = cwiseProduct( *eta_jy_d, *u_zy );
  *ju_zz = cwiseProduct( *eta_jz_d, *u_zz );
 
  thrust::transform( thrust::make_zip_iterator(thrust::make_tuple(ju_xx->begin(), ju_xy->begin(), ju_xz->begin())),
                       thrust::make_zip_iterator(thrust::make_tuple(ju_xx->end(), ju_xy->end(), ju_xz->end())),
                       u_term_stt->begin(),
                       sumOfThree<value_type>() );
  thrust::transform( thrust::make_zip_iterator(thrust::make_tuple(ju_yx->begin(), ju_yy->begin(), ju_yz->begin())),
                       thrust::make_zip_iterator(thrust::make_tuple(ju_yx->end(), ju_yy->end(), ju_yz->end())),
					   u_term_stt->begin() + nx,
                       sumOfThree<value_type>() );
  thrust::transform( thrust::make_zip_iterator(thrust::make_tuple(ju_zx->begin(), ju_zy->begin(), ju_zz->begin())),
                       thrust::make_zip_iterator(thrust::make_tuple(ju_zx->end(), ju_zy->end(), ju_zz->end())),
					   u_term_stt->begin() + 2 * nx,
                       sumOfThree<value_type>() );

  thrust::copy(u_term_stt->begin(), u_term_stt->end(), retMatXd.data());
  return retMatXd;
}

void EffFieldGPU::calcCurlM() {

  *cmx1 = tGradY_cuda->mvp(*mz_d);
  *cmx2 = tGradZ_cuda->mvp(*my_d);
  *cmy1 = tGradZ_cuda->mvp(*mx_d);
  *cmy2 = tGradX_cuda->mvp(*mz_d);
  *cmz1 = tGradX_cuda->mvp(*my_d);
  *cmz2 = tGradY_cuda->mvp(*mx_d);

  thrust::transform( cmx1->begin(), cmx1->end(), cmx2->begin(), curlM->begin(), thrust::minus<value_type>() );
  thrust::transform( cmy1->begin(), cmy1->end(), cmy2->begin(), curlM->begin() + nx, thrust::minus<value_type>() );
  thrust::transform( cmz1->begin(), cmz1->end(), cmz2->begin(), curlM->begin() + 2 * nx, thrust::minus<value_type>() );
}

void EffFieldGPU::setMagDev( MRef& Mag ) {
	copyTimer.start();
	thrust::copy(Mag.col(x).data(), Mag.col(x).data() + nx, mx_d->begin());
	thrust::copy(Mag.col(y).data(), Mag.col(y).data() + nx, my_d->begin());
	thrust::copy(Mag.col(z).data(), Mag.col(z).data() + nx, mz_d->begin());

	thrust::copy( mx_d->begin(), mx_d->end(), mag3_n->begin() );
	thrust::copy( my_d->begin(), my_d->end(), mag3_n->begin() + nx );
	thrust::copy( mz_d->begin(), mz_d->end(), mag3_n->begin() + 2 * nx );
	copyTimer.end();
   copyTimer.add();
}

void EffFieldGPU::setMagDev( const std::vector<value_type>& mag_vec ) {
  nx = mag_vec.size() / 3;
  copyTimer.start();
  thrust::copy( mag_vec.begin(), mag_vec.begin() + nx , mx_d->begin() );
  thrust::copy( mag_vec.begin() + nx, mag_vec.begin() + 2 * nx, my_d->begin() );
  thrust::copy( mag_vec.begin() + 2 * nx, mag_vec.end() , mz_d->begin() );

  thrust::copy( mx_d->begin(), mx_d->end(), mag3_n->begin() );
  thrust::copy( my_d->begin(), my_d->end(), mag3_n->begin() + nx );
  thrust::copy( mz_d->begin(), mz_d->end(), mag3_n->begin() + 2 * nx );
  copyTimer.end();
  copyTimer.add();
}

void EffFieldGPU::setMagDev( const devVecD& mag_vec ) {
  nx = mag_vec.size() / 3;
  copyTimer.start();
  thrust::copy( mag_vec.begin(), mag_vec.begin() + nx , mx_d->begin() );
  thrust::copy( mag_vec.begin() + nx, mag_vec.begin() + 2 * nx, my_d->begin() );
  thrust::copy( mag_vec.begin() + 2 * nx, mag_vec.end() , mz_d->begin() );

  thrust::copy( mx_d->begin(), mx_d->end(), mag3_n->begin() );
  thrust::copy( my_d->begin(), my_d->end(), mag3_n->begin() + nx );
  thrust::copy( mz_d->begin(), mz_d->end(), mag3_n->begin() + 2 * nx );
  copyTimer.end();
  copyTimer.add();
}

void EffFieldGPU::setGradientMatsOnDev( const SpMat& tGradX, const SpMat& tGradY, const SpMat& tGradZ ) {
	tGradX_cuda = std::make_shared<SpMatCUDA>( tGradX );
	tGradY_cuda = std::make_shared<SpMatCUDA>( tGradY );
	tGradZ_cuda = std::make_shared<SpMatCUDA>( tGradZ );
}

void EffFieldGPU::setExchangeMatOnDev(const SpMat& XC_h) {
	xc_cuda = std::make_shared<SpMatCUDA>( XC_h );
	assert (nx == XC_h.rows()) ;
}

void EffFieldGPU::init(int nx) {
	Hxc_unrolled.resize(3*nx);

	mx_d = std::make_shared<dev_vec>(nx);
	my_d = std::make_shared<dev_vec>(nx);
	mz_d = std::make_shared<dev_vec>(nx);
	heffx_d = std::make_shared<dev_vec>(nx);
	heffy_d = std::make_shared<dev_vec>(nx);
	heffz_d = std::make_shared<dev_vec>(nx);
	dxdt = std::make_shared<dev_vec>(3 * nx);
	mag3_n = std::make_shared<dev_vec>(3 * nx);

	dmi3 = std::make_shared<dev_vec>(3 * nx);
	cmx1 = std::make_shared<dev_vec>(nx);
	cmx2 = std::make_shared<dev_vec>(nx);
	cmy1 = std::make_shared<dev_vec>(nx);
	cmy2 = std::make_shared<dev_vec>(nx);
	cmz1 = std::make_shared<dev_vec>(nx);
	cmz2 = std::make_shared<dev_vec>(nx);
	curlM = std::make_shared<dev_vec>(3 * nx);
	surfTerm = std::make_shared<dev_vec>(3 * nx);
	Hxcx_d = std::make_shared<dev_vec>(nx);
	Hxcy_d = std::make_shared<dev_vec>(nx);
	Hxcz_d = std::make_shared<dev_vec>(nx);
	hani_x = std::make_shared<dev_vec>(nx);
	hani_y = std::make_shared<dev_vec>(nx);
	hani_z = std::make_shared<dev_vec>(nx);
	kx_mx = std::make_shared<dev_vec>(nx);
	ky_my = std::make_shared<dev_vec>(nx);
	kz_mz = std::make_shared<dev_vec>(nx);
	tmp0 = std::make_shared<dev_vec>(nx);
	tmp1 = std::make_shared<dev_vec>(nx);
	retVecLLG.resize(3*nx);
//	retVecLLG_d.resize(3*nx);
	retMatXd.resize(nx,3);
}

void EffFieldGPU::setSTTDataOnDevice(const SpMat& GradX, const SpMat& GradY, const SpMat& GradZ,
				     const VectorXd& eta_jx, const VectorXd& eta_jy, const VectorXd& eta_jz) {

	GradX_cuda = std::make_shared<SpMatCUDA>( GradX );
	GradY_cuda = std::make_shared<SpMatCUDA>( GradY );
	GradZ_cuda = std::make_shared<SpMatCUDA>( GradZ );
	ustt_d = std::make_shared<dev_vec>( 3 * nx );
	eta_jx_d = std::make_shared<dev_vec>(nx);
	eta_jy_d = std::make_shared<dev_vec>(nx);
	eta_jz_d = std::make_shared<dev_vec>(nx);
	thrust::copy( eta_jx.data(), eta_jx.data() + nx, eta_jx_d->begin() );
	thrust::copy( eta_jy.data(), eta_jy.data() + nx, eta_jy_d->begin() );
	thrust::copy( eta_jz.data(), eta_jz.data() + nx, eta_jz_d->begin() );
	u_xx = std::make_shared<dev_vec>(nx);
	u_xy = std::make_shared<dev_vec>(nx);
	u_xz = std::make_shared<dev_vec>(nx);
	u_yx = std::make_shared<dev_vec>(nx);
	u_yy = std::make_shared<dev_vec>(nx);
	u_yz = std::make_shared<dev_vec>(nx);
	u_zx = std::make_shared<dev_vec>(nx);
	u_zy = std::make_shared<dev_vec>(nx);
	u_zz = std::make_shared<dev_vec>(nx);

	ju_xx = std::make_shared < dev_vec >(nx);
	ju_xy = std::make_shared < dev_vec >(nx);
	ju_xz = std::make_shared < dev_vec >(nx);
	ju_yx = std::make_shared < dev_vec >(nx);
	ju_yy = std::make_shared < dev_vec >(nx);
	ju_yz = std::make_shared < dev_vec >(nx);
	ju_zx = std::make_shared < dev_vec >(nx);
	ju_zy = std::make_shared < dev_vec >(nx);
	ju_zz = std::make_shared < dev_vec >(nx);

	u_term_stt = std::make_shared < dev_vec >( 3 * nx );
}

void EffFieldGPU::setDMIdata(const VectorXd& D, const VectorXd& invNodeVol,const MatrixXd& nv_n, const VectorXd& nodesurface_h) {
	assert(nx == D.size());
	dev_vec nv_x(nx), nv_y(nx), nv_z(nx);
	dev_vec dmival(nx), invnodevol(nx);
	dev_vec nodesurface(nx);
	thrust::copy(nodesurface_h.data(), nodesurface_h.data() + nx, nodesurface.begin());
	thrust::copy(nv_n.col(x).data(), nv_n.col(x).data() + nx, nv_x.begin());
	thrust::copy(nv_n.col(y).data(), nv_n.col(y).data() + nx, nv_y.begin());
	thrust::copy(nv_n.col(z).data(), nv_n.col(z).data() + nx, nv_z.begin());
	nv_surf_x = std::make_shared<dev_vec>(nx);
	nv_surf_y = std::make_shared<dev_vec>(nx);
	nv_surf_z = std::make_shared<dev_vec>(nx);
	*nv_surf_x = cwiseProduct(nv_x, nodesurface);
	*nv_surf_y = cwiseProduct(nv_y, nodesurface);
	*nv_surf_z = cwiseProduct(nv_z, nodesurface);

	thrust::copy(invNodeVol.data(), invNodeVol.data() + nx, invnodevol.begin());
	thrust::copy(D.data(), D.data() + nx, dmival.begin());
	dev_vec dmi_factor = cwiseProduct(invnodevol, dmival);
	dev_vec dmi_fac3 = dmi_factor;
	dmi_fac3.reserve(2 * nx);
	dmi_fac3.insert(dmi_fac3.end(), dmi_factor.begin(), dmi_factor.end());
	dmi_fac3.insert(dmi_fac3.end(), dmi_factor.begin(), dmi_factor.end()); // yes, do it twice.
	dmi_fac = std::make_shared<dev_vec>(dmi_fac3.size());
	*dmi_fac = dmi_fac3;
};

struct LLG_functor {
  const value_type alpha;
  LLG_functor( const value_type alpha_ )  : alpha(alpha_) {}
    template < class Tuple >
    __device__    
    void operator() ( Tuple t ) const  {

      Matrix< value_type , 3 , 1> M, Heff, MxH, dMdt;
      M(x) = thrust::get<0>( t );
      M(y) = thrust::get<1>( t );
      M(z) = thrust::get<2>( t );

      Heff(x) = thrust::get<6>(t);
      Heff(y) = thrust::get<7>(t);
      Heff(z) = thrust::get<8>(t);

      MxH = M.cross(Heff);

      dMdt = -MxH - alpha * M.cross(MxH);
      thrust::get<3>(t) = dMdt(x);
      thrust::get<4>(t) = dMdt(y);
      thrust::get<5>(t) = dMdt(z);
    }
};

struct LLGnoPrec_functor {
  const value_type alpha;
  LLGnoPrec_functor( const value_type alpha_ )  : alpha(alpha_) {}
    template < class Tuple >
    __device__    
    void operator() ( Tuple t ) const  {
      Matrix< value_type , 3 , 1> M, Heff, MxH, dMdt;

      M(x) = thrust::get<0>( t );
      M(y) = thrust::get<1>( t );
      M(z) = thrust::get<2>( t );

      Heff(x) = thrust::get<6>( t );
      Heff(y) = thrust::get<7>( t );
      Heff(z) = thrust::get<8>( t );

      dMdt = - alpha * M.cross(M.cross(Heff));
      thrust::get<3>(t) = dMdt(x);
      thrust::get<4>(t) = dMdt(y);
      thrust::get<5>(t) = dMdt(z);
    }
};

std::vector<value_type> EffFieldGPU::LLG_noPrec_hst(MRef &Heff,	const value_type alpha) {
	dxdt = LLG_noPrec_dev(Heff, alpha);
	copyTimer.start();
	thrust::copy(dxdt->begin(), dxdt->end(), retVecLLG.begin());
	copyTimer.add();
	return retVecLLG;
}

std::shared_ptr<dev_vec> EffFieldGPU::LLG_noPrec_dev(MRef &Heff, const value_type alpha) {
  thrust::copy( Heff.col(x).data(), Heff.col(x).data() + nx, heffx_d->begin() );
  thrust::copy( Heff.col(y).data(), Heff.col(y).data() + nx, heffy_d->begin() );
  thrust::copy( Heff.col(z).data(), Heff.col(z).data() + nx, heffz_d->begin() );

  thrust::for_each(
		  thrust::make_zip_iterator(
				  thrust::make_tuple(
						  mx_d->begin(),  my_d->begin(),  mz_d->begin() ,
						  dxdt->begin(),
						  dxdt->begin() + nx,
						  dxdt->begin() + 2 * nx,
						  heffx_d->begin(), heffy_d->begin(), heffz_d->begin() )) ,
		  thrust::make_zip_iterator(
				  thrust::make_tuple(
						  mx_d->end(),  my_d->end(), mz_d->end(),
						  dxdt->begin() + nx ,
						  dxdt->begin() + 2 * nx,
						  dxdt->begin() + 3 * nx,
						  heffx_d->end(),  heffy_d->end(), heffz_d->end() )),
				  LLGnoPrec_functor(alpha) );
   return dxdt;
}



std::vector<value_type> EffFieldGPU::ClassicLLG_hst(MRef& Heff, const value_type alpha) {
	dxdt = ClassicLLG_dev(Heff, alpha);
	copyTimer.start();
	thrust::copy(dxdt->begin(), dxdt->end(), retVecLLG.begin());
	copyTimer.add();
	return retVecLLG;
}


std::shared_ptr<dev_vec> EffFieldGPU::ClassicLLG_dev(MRef& Heff, const value_type alpha) {
	copyTimer.start();
	thrust::copy(Heff.col(x).data(), Heff.col(x).data() + nx, heffx_d->begin());
	thrust::copy(Heff.col(y).data(), Heff.col(y).data() + nx, heffy_d->begin());
	thrust::copy(Heff.col(z).data(), Heff.col(z).data() + nx, heffz_d->begin());
	copyTimer.end();
	copyTimer.add();

	thrust::for_each(
			thrust::make_zip_iterator(
					thrust::make_tuple(mx_d->begin(), my_d->begin(), mz_d->begin(),
							dxdt->begin(),
							dxdt->begin() + nx,
							dxdt->begin() + 2 * nx,
							heffx_d->begin(), heffy_d->begin(), heffz_d->begin())),
			thrust::make_zip_iterator(
					thrust::make_tuple(mx_d->end(), my_d->end(), mz_d->end(),
							dxdt->begin() + nx,
							dxdt->begin() + 2 * nx,
							dxdt->begin() + 3 * nx, heffx_d->end(),
							heffy_d->end(), heffz_d->end())),
			LLG_functor(alpha));
	return dxdt;
}



struct STT_functor {
  const value_type alpha;
  const value_type beta;
  STT_functor( const value_type alpha_, const value_type beta_ )  : alpha(alpha_), beta(beta_) {}
    template < class Tuple >
    __device__    
    void operator() ( Tuple t ) const  {

      Matrix< value_type , 3 , 1> M, Ustt, MxU, dMdt;
      M(x) = thrust::get<0>( t );
      M(y) = thrust::get<1>( t );
      M(z) = thrust::get<2>( t );

      Ustt(x) = thrust::get<6>( t );
      Ustt(y) = thrust::get<7>( t );
      Ustt(z) = thrust::get<8>( t );

      MxU = M.cross(Ustt);

      dMdt = -(beta - alpha) * MxU - (1 + alpha * beta) * M.cross(MxU);
      thrust::get<3>( t ) = dMdt(x);
      thrust::get<4>( t ) = dMdt(y);
      thrust::get<5>( t ) = dMdt(z);
    }
};


std::shared_ptr<dev_vec> EffFieldGPU::STT_term_LLG_dev(MRef &Ustt, const value_type alpha, const value_type beta) {

	thrust::copy(Ustt.col(x).data(), Ustt.col(z).data() + nx, ustt_d->begin());
	thrust::for_each(
			thrust::make_zip_iterator(
					thrust::make_tuple(mx_d->begin(), my_d->begin(), mz_d->begin(),
							dxdt->begin(), dxdt->begin() + nx, dxdt->begin() + 2 * nx,
							ustt_d->begin(), ustt_d->begin() + nx, ustt_d->begin() + 2 * nx)),
			thrust::make_zip_iterator(
					thrust::make_tuple(mx_d->end(), my_d->end(), mz_d->end(),
							dxdt->begin() + nx, dxdt->begin() + 2 * nx,	dxdt->begin() + 3 * nx,
							ustt_d->begin() + nx, ustt_d->begin() + 2 * nx, ustt_d->begin() + 3 * nx)),
			STT_functor(alpha, beta));
	return dxdt;
}


std::vector<value_type> EffFieldGPU::STT_term_LLG_hst(MRef &Ustt, const value_type alpha, const value_type beta) {
	dxdt = STT_term_LLG_dev(Ustt, alpha, beta);
	thrust::copy(dxdt->begin(), dxdt->end(), retVecLLG.begin());
	return retVecLLG;
}

struct torque {
  __host__ __device__
  value_type operator() (const thrust::tuple<value_type, value_type, value_type, value_type, value_type, value_type> T ) const  {

    Matrix<value_type , 3 , 1> M, H;
    thrust::tie(M(x), M(y), M(z), H(x), H(y), H(z)) = T;
    return (M.cross(H)).norm();
  }
};

value_type EffFieldGPU::MaxTorque(MRef& Heff) {

	copyTimer.start();
	thrust::copy( Heff.col(x).data(), Heff.col(x).data() + nx, heffx_d->begin() );
	thrust::copy( Heff.col(y).data(), Heff.col(y).data() + nx, heffy_d->begin() );
	thrust::copy( Heff.col(z).data(), Heff.col(z).data() + nx, heffz_d->begin() );
	copyTimer.end();
	copyTimer.add();
	value_type mt = thrust::transform_reduce(
				thrust::make_zip_iterator( thrust::make_tuple(
						mx_d->begin(), my_d->begin(), mz_d->begin(),
						heffx_d->begin(), heffy_d->begin(), heffz_d->begin())),
				thrust::make_zip_iterator( thrust::make_tuple(
						mx_d->end(), my_d->end(), mz_d->end(),
						heffx_d->end(), heffy_d->end(), heffz_d->end())),
						torque(),
						0.0, thrust::maximum<value_type>());
	return mt;
}

Matrix<double, Dynamic,Dynamic> EffFieldGPU::CubicAnisotropyField()	{
	MatrixXd retHcub = MatrixXd::Zero(nx, 3);
	// to be completed

	dev_vec dea_x(nx), dea_y(nx), dea_z(nx);
	dev_vec alf_x(nx), alf_y(nx), alf_z(nx);
	dev_vec a0(nx), a1(nx), a2(nx);
/*
    Vector3d dea;
	for (int i = 0; i < nx; ++i) {
		Vector3d alf = cubicAxes[i].transpose() * Mag.row(i).transpose();
		for (int j = 0; j < 3; ++j) {
			double a0 = alf(j % 3);
			double a1 = alf((j + 1) % 3);
			double a2 = alf((j + 2) % 3);
		    dea(j) = -2. * a0 * (Kc1(i) * (a1 * a1 + a2 * a2) + Kc2(i) * a1 * a1 * a2 * a2);
		}
		Hcub.row(i) = cubicAxes[i] * dea;
	}
*/
	return retHcub;
}


Matrix<double, Dynamic, Dynamic> EffFieldGPU::UniaxialAnisotropyField() {

	*kx_mx = cwiseProduct( *kuAxis_x, *mx_d );
	*ky_my = cwiseProduct( *kuAxis_y, *my_d );
	*kz_mz = cwiseProduct( *kuAxis_z, *mz_d );

	thrust::transform(
			thrust::make_zip_iterator(thrust::make_tuple(kx_mx->begin(), ky_my->begin(), kz_mz->begin())),
	                       thrust::make_zip_iterator(thrust::make_tuple(kx_mx->end(), ky_my->end(), kz_mz->end())),
	                       tmp0->begin(),
	                       sumOfThree<value_type>() );
	*tmp1 = cwiseProduct(*kuVal,*tmp0);
	thrust::transform(tmp1->begin(), tmp1->end(), tmp1->begin(), 2. * thrust::placeholders::_1);

	*hani_x = cwiseProduct( *kuAxis_x, *tmp1 );
	*hani_y = cwiseProduct( *kuAxis_y, *tmp1 );
	*hani_z = cwiseProduct( *kuAxis_z, *tmp1 );

	thrust::copy( hani_x->begin(), hani_x->end(), retMatXd.col(x).data() );
	thrust::copy( hani_y->begin(), hani_y->end(), retMatXd.col(y).data() );
	thrust::copy( hani_z->begin(), hani_z->end(), retMatXd.col(z).data() );
	return retMatXd;
}


struct Cross_GPU {
    template < class Tuple >
    __device__
    void operator() ( Tuple t ) const  {
      Matrix<value_type , 3 , 1> a, b, axb;

      a(x) = thrust::get<0>( t );
      a(y) = thrust::get<1>( t );
      a(z) = thrust::get<2>( t );

      b(x) = thrust::get<3>( t );
      b(y) = thrust::get<4>( t );
      b(z) = thrust::get<5>( t );

      axb = a.cross(b);
      thrust::get<6>( t ) = axb(x);
      thrust::get<7>( t ) = axb(y);
      thrust::get<8>( t ) = axb(z);
    }
};


struct DMI_OP {
	template<class Tuple>
	__host__ __device__
	value_type operator()(Tuple t) const {
		value_type v1, v2, v3;
		thrust::tie(v1, v2, v3) = t;
		return (v1 - 2. * v2) * v3;
	}
};


Matrix<double, Dynamic, Dynamic> EffFieldGPU::DMIField() {
//	Hdmi = (surfIntDMI(Mag) - 2. * CurlM(Mag)).array().colwise() * invNodeVol.cwiseProduct(D).array(); //

	calcCurlM();

	thrust::for_each(
			thrust::make_zip_iterator(
					thrust::make_tuple(mx_d->begin(), my_d->begin(), mz_d->begin(),
							nv_surf_x->begin(), nv_surf_y->begin(), nv_surf_z->begin(),
							surfTerm->begin(), surfTerm->begin() + nx, surfTerm->begin() + 2 * nx)),
			thrust::make_zip_iterator(
					thrust::make_tuple(mx_d->end(), my_d->end(), mz_d->end(),
							nv_surf_x->end(), nv_surf_y->end(), nv_surf_z->end(),
							surfTerm->begin() + nx, surfTerm->begin() + 2 * nx, surfTerm->begin() + 3 * nx) ),
							Cross_GPU());

	thrust::transform(
			thrust::make_zip_iterator(
					thrust::make_tuple(surfTerm->begin(), curlM->begin(),dmi_fac->begin())),
			thrust::make_zip_iterator(
					thrust::make_tuple(surfTerm->end(), curlM->end(), dmi_fac->end())),
					dmi3->begin(),
					DMI_OP());

	thrust::host_vector<value_type> ret_h = *dmi3;
	return Map<Matrix<value_type, Dynamic, Dynamic> >(ret_h.data(), nx, 3);
}


struct normalize {
	template<class Tuple>
	__host__ __device__
	void operator()( Tuple t) {
		Matrix<value_type , 3 , 1> M;
		thrust::tie( M(x), M(y), M(z) ) = t;
		M.normalize();
		thrust::get<x>(t) = M(x);
		thrust::get<y>(t) = M(y);
		thrust::get<z>(t) = M(z);
	}
};



void EffFieldGPU::NormalizeMag( MatrixXd& Mag, int nx ) {
	copyTimer.start();
	thrust::copy(Mag.col(x).data(), Mag.col(x).data() + nx, mx_d->begin());
	thrust::copy(Mag.col(y).data(), Mag.col(y).data() + nx, my_d->begin());
	thrust::copy(Mag.col(z).data(), Mag.col(z).data() + nx, mz_d->begin());
	copyTimer.end();
	copyTimer.add();

	thrust::for_each(
				thrust::make_zip_iterator(
						thrust::make_tuple(mx_d->begin(), my_d->begin(), mz_d->begin())),
				thrust::make_zip_iterator(
						thrust::make_tuple(mx_d->end(), my_d->end(), mz_d->end())),
				normalize());
	copyTimer.start();
	thrust::copy(mx_d->begin(), mx_d->end(), Mag.col(x).data());
	thrust::copy(my_d->begin(), my_d->end(), Mag.col(y).data());
	thrust::copy(mz_d->begin(), mz_d->end(), Mag.col(z).data());
	thrust::copy(mx_d->begin(), mx_d->end(), mag3_n->begin());
	thrust::copy(my_d->begin(), my_d->end(), mag3_n->begin() + nx);
	thrust::copy(mz_d->begin(), mz_d->end(), mag3_n->begin() + 2 * nx);
	copyTimer.end();
	copyTimer.add();
}

void EffFieldGPU::displayTimer() {
	  std::cout << "GPU time for copying data [s]:\t" << copyTimer.durationInMus() / 1.e6 << std::endl;
}

