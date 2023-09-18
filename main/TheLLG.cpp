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
 * TheLLG.cpp
 *
 *  Created on: May 10, 2017
 *      Author: riccardo
 */
///////////////

#include "TheLLG.h"
#include "auxiliaries.h"
#include "PhysicalConstants.h"
#include <iostream>
#include <Eigen/Dense>
#include <cmath>
#include <algorithm>
#include <boost/numeric/odeint.hpp>

#ifdef _OPENMP 
#include <omp.h>
#include <boost/numeric/odeint/external/openmp/openmp.hpp>
#endif


using namespace Eigen;
using namespace boost::numeric::odeint;

static double odeTime;
enum Coords {x, y, z};

// constructor
TheLLG::TheLLG(SimulationData& sd, const MeshData& msh, int LLGVersion) :
		nx(sd.nx),
		NodeVolume(msh.NodeVolume), alpha(sd.alpha), Js(sd.Js),
		hp(sd.Hp), hl(sd.Hl),
		freezeDemag(true)
{
	useGPU = sd.useGPU;
#ifdef USE_CUDA 
	if (useGPU) {
		gpucalc = std::make_shared<EffFieldGPU>(nx);
	}
#endif
	efc = std::make_shared<EffFieldCalc>(sd, msh);
	odeTime = 0.;
	refresh_time = 0.1; // ToDo: consider replacing hard-coded value with user input dat 
	copyTimer.reset();
	heffTimer.reset();
	hl.setValueFunction();
	selectLLGType(LLGVersion);
	invNodeVol = NodeVolume.cwiseInverse();
	Hku   = MatrixXd::Zero(nx, 3);
	Hsurf = MatrixXd::Zero(nx, 3);
	Hexc  = MatrixXd::Zero(nx, 3);
	Hcub  = MatrixXd::Zero(nx, 3);
	Heff  = MatrixXd::Zero(nx, 3);
	Hdmi  = MatrixXd::Zero(nx, 3);
	Hpls  = MatrixXd::Zero(nx, 3);
	Hswp  = MatrixXd::Zero(nx, 3);
	Hdem  = MatrixXd::Zero(nx, 3);
	Hloc  = MatrixXd::Zero(nx, 3);
	ret   = MatrixXd::Zero(nx, 3);
	ret_vec.resize(ret.size());
	useDMI = !sd.D.isZero();
	if (useGPU) {
#ifdef USE_CUDA 
		SpMat XC_field_OP = (-2.) * msh.stiff;
		XC_field_OP = (sd.A.cwiseProduct(invNodeVol)).asDiagonal() * XC_field_OP;
		std::cout << "setting exchange matrix on device" << std::endl;
		gpucalc->setExchangeMatOnDev(XC_field_OP);
		if (useDMI) {
			std::cout << "setting gradient matrices on device" << std::endl;
			gpucalc->setGradientMatsOnDev(-msh.tGradX, -msh.tGradY, -msh.tGradZ);
			gpucalc->setDMIdata(sd.D, invNodeVol, msh.nv_nx, msh.nodeArea);
		}
#endif
	}
	useCubic =  !( sd.Kc1.isZero() && sd.Kc2.isZero() );
	useUniaxial = !sd.Ku1.isZero();
	useSurfaceAnisotropy = !sd.Ks.isZero();
	if (useUniaxial && useGPU) {
	#ifdef USE_CUDA
		gpucalc->setUniaxialAnisotropy(sd.Kuni, sd.Ku1);
	#endif
	}
	realTimeScale = sd.psTimeScaleFactor();
	invJs = Js.cwiseInverse().unaryExpr([](double v) {return std::isfinite(v)? v : 0.0;});
	JsVol = Js.cwiseProduct(NodeVolume);
	if (hp.pulseIsUsed || hp.sweepIsUsed) setHdynField();
	totalVolume = NodeVolume.sum();
	selectEffectiveFields();
}

///////////////////////////// OPERATOR () FOR ODEINT AND CVODE /////////////////////////////////

void TheLLG::operator()(const state_type& state, state_type& dxdt, const double theTime /*time*/) {
	timeInPs = theTime / realTimeScale;
	odeTime = timeInPs;
//	std::cout << "In TheLLG::operator() [CPU] : timeInPs = " << timeInPs << std::endl;
	dxdt = selectedLLGType(state);
}

//////////////////////////// LLG TYPE SELECTION ///////////////////////////////////

void TheLLG::selectLLGType(int choice) {
	enum choices { noPrec , usualLLG , STT };
	useSTT = false;
	if (useGPU) {
#ifdef USE_CUDA
		selectLLGTypeGPU(choice);
#endif
	}
	if (choice == noPrec) {
		selectedLLGType = [this](const state_type& m)->state_type { return noPrecession(m); };
	} else if (choice == STT ) {
		useSTT = true;
		selectedLLGType = [this](const state_type& m)->state_type { return sttDynamics(m); };
	} else {
		selectedLLGType = [this](const state_type& m)->state_type { return classicVersion(m); };
	}
}


/////////////////////////// LLG VERSIONS //////////////////////////


state_type TheLLG::classicVersion(const state_type &mag_vec) {
	if (useGPU) {
	#ifdef USE_CUDA
			gpucalc->setMagDev(mag_vec);
	#endif
		}  // requires ColMajor storage for Mag
	Map<const MatrixXd_CM> Mag(mag_vec.data(), nx, 3);
	evaluateAllEffectiveFields(Mag);
	Heff = totalEffectiveField();
	if (useGPU) {
#ifdef USE_CUDA
// This is invoked in the option "integrate LLG on CPU while using GPU for eff. field calculation"
		ret_vec = gpucalc->ClassicLLG_hst(Heff, alpha);
#endif
	} else {
		const MatrixXd MxH = cross(Mag, Heff);
//		ret = -MxH - alpha * cross(Mag, MxH);
		/*
		 for (int i = 0; i < nx; ++i) {
		 Vector3d M = Mag.row(i);
		 Vector3d Hloc = Heff.row(i);
		 ret.row(i) = - M.cross(Hloc) - alpha * M.cross(M.cross(Hloc));
		 }
		 */
//		Map<MatrixXd_CM>(ret_vec.data(), ret.rows(), 3) = ret;
		Map<MatrixXd_CM>(ret_vec.data(), ret.rows(), 3) = -MxH - alpha * cross(Mag, MxH);
	}
	return ( ret_vec );
}

state_type TheLLG::noPrecession(const state_type &mag_vec) {
	if (useGPU) {
#ifdef USE_CUDA
		gpucalc->setMagDev(mag_vec);
#endif
	}
	Map<const MatrixXd_CM> Mag(mag_vec.data(), nx, 3);
	evaluateAllEffectiveFields(Mag);
	Heff = totalEffectiveField();
	if (useGPU) {
#ifdef USE_CUDA
		ret_vec = gpucalc->LLG_noPrec_hst(Heff, alpha);
#endif
	} else {
//		ret = -alpha * cross(Mag, cross(Mag, Heff)); // host version
//		Map<MatrixXd_CM>(ret_vec.data(), ret.rows(), 3) = ret;
		Map<MatrixXd_CM>(ret_vec.data(), ret.rows(), 3) = -alpha * cross(Mag, cross(Mag, Heff));
	}
	return (ret_vec);
}

state_type TheLLG::sttDynamics(const state_type &mag_vec) {
	const std::vector<double> LLGpart = classicVersion(mag_vec);
	Map<const MatrixXd_CM> Mag(mag_vec.data(), nx, 3);
	calcUtermSTT(Mag);
	if (useGPU) {
#ifdef USE_CUDA
		ret_vec = gpucalc->STT_term_LLG_hst(stt.Ustt, alpha, stt.beta);
#endif
	} else {
		const MatrixXd MxU = cross(Mag, stt.Ustt);
		ret = -(stt.beta - alpha) * MxU	- (1. + alpha * stt.beta) * cross(Mag, MxU);
		Map<MatrixXd_CM>(ret_vec.data(), ret.rows(), 3) = ret;
	}
	if (stt.pulseIsUsed) {
		double pulseVal = stt.gaussPulseValue(timeInPs);
		std::transform(ret_vec.begin(), ret_vec.end(), ret_vec.begin(), [pulseVal](double& c){ return c * pulseVal;});
	}
	std::transform(ret_vec.begin(), ret_vec.end(), LLGpart.begin(),	ret_vec.begin(), std::plus<double>());
	return (ret_vec);
}


////////////////////// EFFECTIVE FIELDS ////////////////////////////

void TheLLG::selectEffectiveFields() {
	calcEffectiveField.clear();
	calcEffectiveField.emplace_back( [this](MRef &m) { calcExchangeField(m); } );
	if (useUniaxial) calcEffectiveField.emplace_back( [this](MRef &m) { calcUniaxialAnisotropyField(m); } );
	if (useCubic) calcEffectiveField.emplace_back( [this](MRef &m) { calcCubicAnisotropyField(m); } );
	if (useDMI) calcEffectiveField.emplace_back( [this](MRef &m) { calcDMIField(m); } );
	if (useSurfaceAnisotropy) calcEffectiveField.emplace_back( [this](MRef &m) { calcSurfaceAnisotropyField(m); } );
	if (hp.pulseIsUsed) calcEffectiveField.emplace_back( [this](MRef &m){ calcPulseField(m); } );
	if (hp.sweepIsUsed) calcEffectiveField.emplace_back( [this](MRef &m) { calcSweepField(m); } );
	if (hl.isUsed) calcEffectiveField.emplace_back( [this](MRef &m) { calcLocalField(m); } );
}


MatrixXd TheLLG::totalEffectiveField() {
	MatrixXd Htot(nx, 3);
#ifdef _OPENMP
		int i,j;
#pragma omp parallel for private(i) collapse(2)
		for (i = 0; i < nx; ++i) {
			for (j = 0; j < 3; ++j ) {
				Htot(i, j) = (Hexc(i, j) + Hku(i, j) + Hcub(i, j) + Hsurf(i, j) + Hdmi(i, j)) * invJs(i)
						+ Hext(i, j) + Hdem(i,j) + Hpls(i, j) + Hswp(i, j) + Hloc(i,j);
			}
		}
#else
	MatrixXd Heff_tot = (Hexc + Hku + Hcub + Hsurf + Hdmi).array().colwise()* invJs.array();
	Htot = Heff_tot + Hext + Hdem + Hpls + Hswp + Hloc;
#endif
	return Htot;
}


void TheLLG::evaluateAllEffectiveFields( MRef &Mag ) {
	int tasks = calcEffectiveField.size();
	for (int i = 0; i < tasks; ++i) {
		calcEffectiveField[i](Mag);
	}
}


void TheLLG::calcExchangeField(MRef &Mag) {
	heffTimer.start();
	if (useGPU) {
#ifdef USE_CUDA 
		Hexc = gpucalc->ExchangeFieldGPU();
#endif
	} else {
		Hexc = efc->exchangeField(Mag);
	}
	heffTimer.add();
}


void TheLLG::calcDemagField(MRef &Mag) {
	if (odeTime >= next_refresh) {
		next_refresh += refresh_time;
		Hdem = demag->calcField(Mag);
	}
}


void TheLLG::calcUniaxialAnisotropyField(MRef &Mag) {
	if (useGPU) {
#ifdef USE_CUDA
		Hku = gpucalc->UniaxialAnisotropyField();
#endif
	} else {
		Hku = efc->uniaxialAnisotropyField(Mag);
	}
}


void TheLLG::calcCubicAnisotropyField( MRef &Mag) {
	if (useCubic) {
		Hcub = efc->cubicAnisotropyField(Mag);
	}
}


void TheLLG::calcSurfaceAnisotropyField(MRef &Mag) {
	Hsurf = efc->surfaceAnisotropyField(Mag);
}


void TheLLG::calcDMIField(MRef &Mag) {
	if (useGPU) {
#ifdef USE_CUDA
		Hdmi = gpucalc->DMIField();
#endif
	} else {
		Hdmi = efc->dmiField(Mag);
	}
}


void TheLLG::calcPulseField( MRef& ){
	Hpls = hp.pulseH * hp.gaussPulseValue(timeInPs);
}


void TheLLG::calcSweepField( MRef& ){
	Hswp = hp.sweepH * hp.sweepFieldValue(timeInPs);
}


void TheLLG::calcLocalField( MRef& ){
	Hloc = hl.localField(timeInPs);
}

///////////////////////  ENERGIES  ///////////////////////

double TheLLG::getExchangeEnergy( MRef &Mag ) {
	return efc->calcExchangeEnergy(Mag);
}
double TheLLG::getDirectExch(MRef& Mag)  {
	return efc->exchEnergy_direct(Mag);
}



double TheLLG::getZeemanEnergy(MRef &Mag) {
	double zeeEnergy =  -(JsVol.transpose() * (Hext + Hpls + Hswp + Hloc).cwiseProduct(Mag)).sum();
	return zeeEnergy;
}


double TheLLG::getDemagEnergy( MRef &Mag ) {
	if(freezeDemag) calcDemagField(Mag);
	return demag->getDemagEnergy(Mag);
}


double TheLLG::getUniaxialAnisotropyEnergy( MRef &Mag ){
	if ( !useUniaxial ) return 0;
	return efc->uniaxialEnergy(Mag);
}


double TheLLG::getCubicAnisotropyEnergy( MRef &Mag ) {
	if (!useCubic) return 0;
	return efc->cubicAnisotropyEnergy(Mag);
}


double TheLLG::getDMIEnergy( MRef &Mag ) {
	if (!useDMI) return 0;
#ifdef USE_CUDA
	calcDMIField(Mag);
	efc->setDMIField(Hdmi);
#endif
	return efc->dmiEnergy(Mag);
}


double TheLLG::getSurfaceAnisotropyEnergy(MRef &Mag) {
	if (!useSurfaceAnisotropy) return 0.;
	return efc->surfaceAnisotropyEnergy(Mag);
}

double TheLLG::getDirectDMI(MRef& Mag)  {
	return efc->dmiEnergy_direct(Mag);
}


/////////////////////////// OTHER /////////////////////////////////////////////


void TheLLG::outputTimer() {
      std::cout << "time for copying data [s]:\t" << std::setprecision(2) << copyTimer.durationInMus() / 1.e6 << std::endl;
#ifdef USE_CUDA
      gpucalc->displayTimer();
#endif
      std::cout << "time for effective field: " << std::setprecision(2) << heffTimer.durationInMus() / 1.e6 << std::endl;
      std::cout << std::setprecision(-1);
}


void TheLLG::calcUtermSTT(MRef &Mag) {
	if (useGPU) {
#ifdef USE_CUDA
		stt.Ustt = gpucalc->UTermSTT_GPU();
#endif
	} else {
/*
		VectorXd Mx = Mag.col(x);
		VectorXd My = Mag.col(y);
		VectorXd Mz = Mag.col(z);
		stt.Ustt.col(x) = stt.eta_jx.cwiseProduct(stt.gradX * Mx)
				+ stt.eta_jy.cwiseProduct(stt.gradY * Mx)
				+ stt.eta_jz.cwiseProduct(stt.gradZ * Mx);
		stt.Ustt.col(y) = stt.eta_jx.cwiseProduct(stt.gradX * My)
				+ stt.eta_jy.cwiseProduct(stt.gradY * My)
				+ stt.eta_jz.cwiseProduct(stt.gradZ * My);
		stt.Ustt.col(z) = stt.eta_jx.cwiseProduct(stt.gradX * Mz)
				+ stt.eta_jy.cwiseProduct(stt.gradY * Mz)
				+ stt.eta_jz.cwiseProduct(stt.gradZ * Mz);
		*/
// #pragma omp parallel for // SparseMVP is already parallelized by Eigen.
		for (int i = 0; i < 3; ++i) {
			stt.Ustt.col(i) = stt.eta_jx.cwiseProduct(stt.gradX * Mag.col(i))
							+ stt.eta_jy.cwiseProduct(stt.gradY * Mag.col(i))
							+ stt.eta_jz.cwiseProduct(stt.gradZ * Mag.col(i));

		}
	}
}


/////////////////// GETTER  //////////////////////////

Vector3d TheLLG::getMeanH() {
	return  NodeVolume.transpose() * (Hext + Hpls + Hswp + Hloc)  / totalVolume;
}


double TheLLG::getMaxTorque( MRef &Mag ) {
	evaluateAllEffectiveFields( Mag );
	Heff = totalEffectiveField();
	double m_torque = 0;
	if (useGPU) {
#ifdef USE_CUDA
		m_torque = gpucalc->MaxTorque(Heff);
#endif
	} else {
		m_torque = cross(Mag, Heff).rowwise().norm().maxCoeff();
	}
// Note: The max.-torque criterion only includes effective field contributions (no STT torques).
	return m_torque * PhysicalConstants::mu0;
}


MatrixXd fieldFromAngles(int nx, double Habs, double theta_h, double phi_h) {
	MatrixXd Hfield = MatrixXd::Ones(nx, 3); // use other distribution for possible extension of inhomogeneous field
	Hfield.col(x) *= std::sin(theta_h) * std::cos(phi_h);
	Hfield.col(y) *= std::sin(theta_h) * std::sin(phi_h);
	Hfield.col(z) *= std::cos(theta_h);
	Hfield *= ( Habs / PhysicalConstants::mu0 );
	return Hfield;
}


/////////////////////////////////// SETTER ///////////////////////////////////////////

void TheLLG::setTime(double time) {
	timeInPs = time;
}


void TheLLG::setDemagData(const DemagField& demag_) {
	demag = std::make_shared<DemagField>(demag_);
	freezeDemag = false;
	calcEffectiveField.emplace_back([this](MRef &m) { calcDemagField(m); } );
}


void TheLLG::setSTTData(const STT &stt_) {
	stt = stt_;
	stt.setEta(invJs);
	if (useGPU) {
#ifdef USE_CUDA
		std::cout << "setting STT data on device" << std::endl;
		gpucalc->setSTTDataOnDevice(stt.gradX, stt.gradY, stt.gradZ, stt.eta_jx, stt.eta_jy, stt.eta_jz);
#endif
	}
}


void TheLLG::setZeemanField(int nx, double H_stat, double theta_H_deg, double phi_H_deg,
		double H_hys, double theta_Hys_deg, double phi_Hys_deg) {
	double deg_to_rad = PhysicalConstants::pi / 180.;
	MatrixXd Hhys = fieldFromAngles(nx, H_hys, theta_Hys_deg * deg_to_rad, phi_Hys_deg * deg_to_rad);
	MatrixXd Hstat = fieldFromAngles(nx, H_stat, theta_H_deg * deg_to_rad, phi_H_deg * deg_to_rad);
	Hext = Hhys + Hstat;
}


void TheLLG::setHdynField () {
	hp.pulseH = MatrixXd::Zero(nx,3);
	hp.sweepH = MatrixXd::Zero(nx,3);
	if (hp.pulseIsUsed) {
		hp.pulseH = fieldFromAngles(nx, 1. , hp.pulseTheta, hp.pulsePhi);
	}
	if (hp.sweepIsUsed) {
		hp.sweepH = fieldFromAngles(nx, 1. , hp.sweepTheta, hp.sweepPhi);
	}
}


void TheLLG::setMag(MRef &mag_) {
	(void) mag_;
	if (useGPU) {
#ifdef USE_CUDA
		gpucalc->setMagDev(mag_);
#endif
	}
}


void TheLLG::setHdem(MRef &Hdem_) {
	Hdem = Hdem_;
}

//////////////////////////   ODEINT INTEGRATOR  //////////////////////////////////////////////

int TheLLG::IntegrateODEINT(std::vector<double>& mag_vec, double ode_start_t, double ode_end_t, double dt ) {
	std::reference_wrapper<TheLLG> LLGRef = std::ref(*this);
#ifdef _OPENMP 
	typedef runge_kutta_dopri5< state_type , double , state_type , double , openmp_range_algebra > dopri54_omp;
	typedef controlled_runge_kutta< dopri54_omp > dopri54_controlled_omp;
#endif	
	int its;
#ifdef _OPENMP 
//			total_its += integrate_adaptive(adb_stepper(), std::ref(LLG), mag_vec, ode_start_t, ode_start_t + ode_end_t, dt);
//			total_its += integrate_adaptive(dopri54_controlled_omp(), std::ref(LLG), mag_vec, ode_start_t, ode_start_t + ode_end_t, dt, myobs);
			its = integrate_adaptive(dopri54_controlled_omp(), LLGRef, mag_vec, ode_start_t, ode_start_t + ode_end_t, dt);
//			its = integrate_adaptive(dopri54_controlled_omp(), this, mag_vec, ode_start_t, ode_start_t + ode_end_t, dt);
//			total_its += integrate_adaptive(fehlberg78_controlled_omp(), std::ref(LLG), mag_vec, ode_start_t, ode_start_t + ode_end_t, dt);
//#endif
#else
  		its = integrate<double>( LLGRef, mag_vec, ode_start_t, ode_start_t + ode_end_t, dt );
#endif
  		//]
  		//[ Bulirsch-Stoer Stepper:
  		//		total_its += integrate_adaptive(bulirschStoerStepper, std::ref(LLG), mag_vec, ode_start_t, ode_start_t + ode_end_t, dt); // relatively low number of iterations, but overall slower than default
  		//]
  		//[ OpenMp integration:
  		//		total_its += integrate_adaptive(bulirschStoer_omp, std::ref(LLG), mag_vec, ode_start_t, ode_start_t + ode_end_t, dt); // does not compile with boost v.1.58.0. Bug in odeint? -> Works with boost 1.69
  		//		int nsteps = 1;
  		//		total_its += integrate_n_steps(bulirschStoer_omp, std::ref(LLG), mag_vec, ode_start_t,  ode_end_t / nsteps   , nsteps);
  		//]
  		/*
  		//[ Default integrator (Dormand-Prince Method)
  		#ifdef _OPENMP
	    //typedef openmp_state<std::vector<double>> omp_state_type;
  		//			total_its += integrate_adaptive(adb_stepper(), std::ref(LLG), mag_vec, ode_start_t, ode_start_t + ode_end_t, dt);
  		//			total_its += integrate_adaptive(dopri54_controlled_omp(), std::ref(LLG), mag_vec, ode_start_t, ode_start_t + ode_end_t, dt, myobs);
  					total_its += integrate_adaptive(dopri54_controlled_omp(), std::ref(LLG), mag_vec, ode_start_t, ode_start_t + ode_end_t, dt);
  		//			total_its += integrate_adaptive(fehlberg78_controlled_omp(), std::ref(LLG), mag_vec, ode_start_t, ode_start_t + ode_end_t, dt);
  		//#endif
  		#else
  		  		total_its += integrate<double>( std::ref(LLG), mag_vec, ode_start_t, ode_start_t + ode_end_t, dt );
  		#endif
  		*/
  		//	typedef runge_kutta_cash_karp54<state_type> error_stepper_type;
  		//	typedef runge_kutta_fehlberg78<state_type, double, state_type, double, openmp_range_algebra> fehlberg78_omp;
  		//    typedef runge_kutta_dopri5< state_type , double , state_type , double , openmp_range_algebra > dopri54_omp;
  		//	typedef controlled_runge_kutta< dopri54_omp > dopri54_controlled_omp;
  		//	typedef controlled_runge_kutta< fehlberg78_omp > fehlberg78_controlled_omp;
  			// for some strange reason the following is not compatible with GPU calculation. It leads to inconsistent results.
  			// typedef adams_bashforth_moulton< 4 , state_type, double, state_type, double, openmp_range_algebra  > adb_stepper;

  		//	bulirsch_stoer< state_type> bulirschStoerStepper( 1E-6 , 0.0 , 0.0 , 0.0);
  		//	bulirsch_stoer< state_type, double, state_type, double, openmp_range_algebra> bulirschStoer_omp;
  		//	bulirsch_stoer< state_type> bulirschStoerStepper;

		//#ifdef _OPENMP
		//	const int chunk_size =   mag.rows() / omp_get_max_threads();
		//#endif
  		//#ifdef _OPENMP
  		//		omp_set_schedule( omp_sched_static, chunk_size );
  		//#endif
return its;
}

//////////////////////////////////// SUNDIALS INTEGRATOR //////////////////////////////////////////////////////

#define NV_Ith_S(v,i) ( NV_DATA_S(v)[i] ) // for CVODE

int TheLLG::integrateSUNDIALS(std::vector<double>& mag_vec, double ode_start_t, double ode_end_t, double dt ) {
	(void) dt;
	long int its_l, its_nl = 0;
	int flag;
	sunindextype N = 3 * nx;
	N_Vector m;
	copyTimer.start();
#ifndef USE_CVODE_5
	SUNContext sunctx; 
    SUNContext_Create(NULL, &sunctx);
# endif	
#ifdef _OPENMP
    int num_threads = omp_get_max_threads();
//	 int num_threads = omp_get_num_threads();	
	#ifdef USE_CVODE_5
		m = N_VMake_OpenMP(N, mag_vec.data(), num_threads); 
	#else			
		m = N_VMake_OpenMP(N, mag_vec.data(), num_threads, sunctx); 
	#endif	
#else	
	#ifdef USE_CVODE_5
		m = N_VMake_Serial(N, mag_vec.data()); 	
	#else
		m = N_VMake_Serial(N, mag_vec.data(), sunctx); 		
	#endif 
#endif
	copyTimer.add();
	void *cvode_mem = NULL;
	#ifdef USE_CVODE_5
		cvode_mem = CVodeCreate(CV_ADAMS);
	#else
		cvode_mem = CVodeCreate(CV_ADAMS, sunctx); 		
	#endif
//	cvode_mem = CVodeCreate(CV_BDF);
	realtype t0 = ode_start_t;
	flag = CVodeInit(cvode_mem, rhs, t0, m);
	if (flag) return(1);

	realtype abstol = 1.e-6;
	realtype reltol = 1.e-6;
	flag = CVodeSStolerances(cvode_mem, reltol, abstol);
	std::shared_ptr<UserData> data = alloc_user_data(nx);
	flag = CVodeSetUserData(cvode_mem, data.get());

	SUNLinearSolver LS;
	// LS = SUNSPGMR(m, 0, 0);	
	#ifdef USE_CVODE_5
		LS = SUNLinSol_SPGMR(m, PREC_NONE, 0); 	
	#else 
		LS = SUNLinSol_SPGMR(m, PREC_NONE, 0, sunctx); 		
	#endif	
	flag = CVodeSetLinearSolver(cvode_mem, LS, NULL);

	realtype tout = ode_start_t + ode_end_t;
	realtype t = ode_start_t;
	flag = CVode(cvode_mem, tout, m, &t, CV_NORMAL);
	CVodeGetNumNonlinSolvIters(cvode_mem, &its_nl);
	CVodeGetNumLinIters(cvode_mem, &its_l);
	int its = static_cast<int>(its_nl + its_l);

	realtype* res_p = N_VGetArrayPointer(m);
	copyTimer.start();
	std::copy(res_p, res_p + N, mag_vec.begin());
	copyTimer.add();

#ifdef _OPENMP 
	N_VDestroy_OpenMP(m);
#else
	N_VDestroy(m);
#endif
	CVodeFree(&cvode_mem);
	SUNLinSolFree(LS);
	#ifndef USE_CVODE_5
		SUNContext_Free(&sunctx); 
	#endif
	return its;
}


std::shared_ptr<UserData> TheLLG::alloc_user_data(const int nx) {
  UserData d;
  d.llg = std::make_shared<TheLLG>(std::reference_wrapper<TheLLG>(*this));
  d.nx = nx;
  d.ret.reserve(3*nx);
  std::shared_ptr<UserData> data = std::make_shared<UserData>(d);
//  d.llg = std::make_shared<TheLLG>(*this); // OK
//  d.llg = std::make_shared<TheLLG>(std::ref(*this)); // is that better in terms of performance?
//  d.llg = std::make_shared<TheLLG>(std::reference_wrapper<TheLLG>(*this)); // should be the same as above
  return data;
}


int TheLLG::rhs(realtype t, N_Vector u, N_Vector u_dot, void *user_data) {
	(void) t;
	UserData* u_data = (UserData*)user_data;
	int nx = u_data->nx;
	realtype *udata;
//  realtype *dudata;
	u_data->llg->copyTimer.start();

#ifdef _OPENMP 
//  dudata = NV_DATA_OMP(u_dot);
  udata = NV_DATA_OMP(u);
#else
	udata = N_VGetArrayPointer(u);
//  dudata = N_VGetArrayPointer(u_dot);
#endif
	std::vector<double> mag_vec(udata, udata + (3 * nx));
	u_data->llg->copyTimer.add();
	u_data->llg->operator()(mag_vec, u_data->ret, t);
// this->operator()( mag_vec, u_data->ret, t ); // 'this' is not available in static member functions
	int i;
	u_data->llg->copyTimer.start();
//#ifdef _OPENMP
	//int num_threads = omp_get_max_threads();
	//dudata = u_data->ret.data();
	//u_dot = N_VMake_OpenMP(3 * nx, dudata, num_threads);
//#endif

#ifdef _OPENMP
  #pragma omp parallel for
#endif

  for (i = 0; i < 3 * nx; ++i ) {
#ifdef _OPENMP 
    NV_Ith_OMP(u_dot,i) = u_data->ret[i];
#else
    NV_Ith_S(u_dot,i) = u_data->ret[i];
//    dudata[i] = u_data->ret[i];
#endif
    u_data->llg->copyTimer.add();
  }
  return(0); // <--- This is required to signal success
}

//////////////// END SUNDIALS //////////////////////////////


