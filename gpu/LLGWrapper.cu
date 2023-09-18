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
 * LLGWrapper.cu
 *
 *  Created on: Sep 23, 2019
 *      Author: riccardo
 */

// Cannot be included in TheLLG class: Separate translation unit needed for nvcc.

#include "LLGWrapper.h"
#include <functional>
#include <thrust/copy.h>
#include <nvector/nvector_cuda.h>
#include <thrust/transform.h>
#include <thrust/functional.h>
#include <thrust/execution_policy.h>

dev_vec mag_vec_tmp;
dev_vec ret_vec_tmp;

// using namespace boost::numeric::odeint;
//typedef runge_kutta_fehlberg78<dev_vec, double, dev_vec, double> fehlberg78_gpu;
//typedef controlled_runge_kutta< fehlberg78_gpu > fehlberg78_controlled_gpu;


void LLGWrapper::init(const int nx) {
	nx3 = 3 * nx;
	mag_vec_tmp.resize(nx3);
	ret_vec_tmp.resize(nx3);
}


dev_vec TheLLG::sttDynamics_GPU(const dev_vec &mag_vec) {
	dev_vec LLGpart_d = classicVersion_GPU(mag_vec);
	stt.Ustt = gpucalc->UTermSTT_GPU();
	dev_vec ret_vec_d = *gpucalc->STT_term_LLG_dev(stt.Ustt, alpha, stt.beta);
	if (stt.pulseIsUsed) {
		double pulseVal = stt.gaussPulseValue(timeInPs);
		thrust::transform(ret_vec_d.begin(), ret_vec_d.end(), ret_vec_d.begin(), pulseVal * thrust::placeholders::_1);
	}
	thrust::transform(thrust::device, ret_vec_d.begin(), ret_vec_d.end(), LLGpart_d.begin(), ret_vec_d.begin(), thrust::plus<double>());
	return(ret_vec_d);
}

Eigen::MatrixXd TheLLG::effFieldsForGPU(const dev_vec &mag_vec) {
// This should be improved. There are too many copy operations.	
// An overlaod of evaluateAllEffectiveFields() could be prepared
// that would take a device vector as argument.

	std::vector<double> mag_vec_h(3 * nx);
	thrust::copy(mag_vec.begin(), mag_vec.end(), mag_vec_h.begin());
	Eigen::Map<const MatrixXd_CM> Mag(mag_vec_h.data(), nx, 3);
	evaluateAllEffectiveFields(Mag);
	return totalEffectiveField();
// One should further check if it is necessary to return MatrixXd 
// in the individual field calculations. 
// It should be possible to leave the data on the device.
}

dev_vec TheLLG::classicVersion_GPU(const dev_vec &mag_vec) {
	gpucalc->setMagDev(mag_vec);
	Heff = effFieldsForGPU(mag_vec);
	return *gpucalc->ClassicLLG_dev(Heff, alpha);
}


dev_vec TheLLG::noPrecession_GPU(const dev_vec &mag_vec) {
	gpucalc->setMagDev(mag_vec);
	Heff = effFieldsForGPU(mag_vec);
	return *gpucalc->LLG_noPrec_dev(Heff, alpha);
}


template<class deviceType>
void delete_vec( deviceType& x_d ) {
	x_d.clear();
	x_d.shrink_to_fit();
	x_d.~device_vector();
}

LLGWrapper::~LLGWrapper(){
	cudaDeviceSynchronize();
	delete_vec(mag_vec_tmp);
	delete_vec(ret_vec_tmp);
}

void TheLLG::operator()(const thrust::device_vector<double>& mag_d, thrust::device_vector<double>& dxdt_d, const double theTime /*time*/) {
	timeInPs = theTime / realTimeScale;
	dxdt_d = selectedLLGType_GPU(mag_d);
}


void TheLLG::selectLLGTypeGPU(int choice) {
	enum choices { noPrec , usualLLG , STT };
	useSTT = false;
		if (choice == noPrec) {
			selectedLLGType_GPU = [this](const dev_vec &m) -> dev_vec { return noPrecession_GPU(m); };
		} else if (choice == STT ) {
			useSTT = true;
			selectedLLGType_GPU = [this](const dev_vec &m) -> dev_vec { return sttDynamics_GPU(m); };
		} else	{
			selectedLLGType_GPU = [this](const dev_vec &m) -> dev_vec {	return classicVersion_GPU(m); };
		}
}


int TheLLG::IntegrateOnGPU(std::vector<double>& mag_vec, double ode_start_t, double ode_end_t, double dt ) {
// Not referenced any more. Used to be called from TheSimulation::start() but was removed due to THRUST / ODEINT incompatibility
	std::reference_wrapper<TheLLG>  LLGRef = std::ref(*this);
	copyTimer.start();
	mag_vec_tmp = mag_vec; // copy to device vector
    copyTimer.add();
//    int its = boost::numeric::odeint::integrate<double>(LLGRef, mag_vec_tmp, ode_start_t, ode_start_t + ode_end_t, dt );
    copyTimer.start();
    thrust::copy(mag_vec_tmp.begin(), mag_vec_tmp.end(), mag_vec.begin());
    copyTimer.add();
//    return its;
	return 0;
}


int TheLLG::gpuODE(std::vector<double>& mag_vec, double ode_start_t, double ode_end_t, double dt ) {
	(void) dt;
	long int its_l, its_nl = 0;
	std::shared_ptr<UserData> data = alloc_user_data(nx);
	int flag;
	sunindextype N = 3 * nx;
    copyTimer.start();
//
//	int* major = (int*)malloc(sizeof(int));
//	int* minor = (int*)malloc(sizeof(int)); 
//	int* patch = (int*)malloc(sizeof(int));
//	int major ;
//	int minor ;
//	int patch ;
//	int len;
//	char label;
//	char* label = (char *)malloc( 20 *sizeof(char));;
//	SUNDIALSGetVersionNumber(major, minor, patch, label, len);
//	int SUNDIALSGetVersionNumber(major, minor, 	patch, 	&label, len);
//	std::cout << "CVODE version : " << major << "." << minor << "." << patch << std::endl;
//	std::exit(1);
//	


	mag_vec_tmp = mag_vec;  // copy to thrust::device_vector
#ifndef OLD_CVODE_VERSION
	SUNContext sunctx; 
    SUNContext_Create(NULL, &sunctx); 
	N_Vector m_gpu = N_VMake_Cuda(N, mag_vec.data(), thrust::raw_pointer_cast(mag_vec_tmp.data()), sunctx); 
#else	
    N_Vector m_gpu = N_VMake_Cuda(N, mag_vec.data(), thrust::raw_pointer_cast(mag_vec_tmp.data())); 
#endif	
    copyTimer.add();

// copy to N_V
//	cudaMemcpy(mdata, thrust::raw_pointer_cast(mag3_d.data()), N * sizeof(realtype), cudaMemcpyDeviceToDevice);
	void *cvode_mem = NULL;	
	#ifndef OLD_CVODE_VERSION
		cvode_mem = CVodeCreate(CV_ADAMS, sunctx); 
	#else
		cvode_mem = CVodeCreate(CV_ADAMS); 
	#endif
//	cvode_mem = CVodeCreate(CV_BDF);
	flag = CVodeInit(cvode_mem, rhs_d, ode_start_t, m_gpu);
	
	// CV_BDF configuration:
	//	realtype abstol = 1.e-8;
	//	realtype reltol = 1.e-5;

	// CV_ADAMS configuration:
	realtype abstol = 1.e-6;
	realtype reltol = 1.e-6;
	flag = CVodeSStolerances(cvode_mem, reltol, abstol);
	flag = CVodeSetUserData(cvode_mem, data.get());

	SUNLinearSolver LS;
	#if OLD_CVODE_VERSION
		LS = SUNLinSol_SPGMR(m_gpu, PREC_NONE, 0); 	
	#else		
		LS = SUNLinSol_SPGMR(m_gpu, PREC_NONE, 0, sunctx); 	
	#endif
	flag = CVodeSetLinearSolver(cvode_mem, LS, NULL);

	realtype tout = ode_end_t + ode_start_t;
	realtype t = ode_start_t;

	flag = CVode(cvode_mem, tout, m_gpu, &t, CV_NORMAL);
	if (flag) {std::cerr << "Warning: GPU integration failed." << std::endl; }
	CVodeGetNumNonlinSolvIters(cvode_mem, &its_nl);
	CVodeGetNumLinIters(cvode_mem, &its_l);
	int its = static_cast<int>(its_nl + its_l); // iterations
	double* res_p = N_VGetDeviceArrayPointer_Cuda(m_gpu);
	copyTimer.start();
	thrust::copy(res_p, res_p + N, mag_vec_tmp.begin()); // copy from N_V
	thrust::copy(mag_vec_tmp.begin(), mag_vec_tmp.end(), mag_vec.begin()); // copy to std::vector -- returned by ref.
	copyTimer.add();
//  The magnetization needs to be transfered to EffFieldGPU in order to calculate max.torque
	gpucalc->setMagDev(mag_vec_tmp);
//	N_VCopyFromDevice_Cuda(m);
//	double * res_p = N_VGetHostArrayPointer_Cuda(m);
//	thrust::copy(res_p, res_p + N, mag_vec.begin());
	N_VDestroy_Cuda(m_gpu);
	CVodeFree(&cvode_mem);
	SUNLinSolFree(LS);
	#ifndef OLD_CVODE_VERSION
		SUNContext_Free(&sunctx); 
	#endif 
	return its;
}


int TheLLG::rhs_d(realtype t, N_Vector u, N_Vector u_dot, void *user_data) {
	(void) t;
	UserData *u_data;
	u_data = (UserData*) user_data;
	sunindextype N = 3 * u_data->nx;
	realtype *dudata = N_VGetDeviceArrayPointer_Cuda(u_dot);
	realtype *udata =  N_VGetDeviceArrayPointer_Cuda(u);
	u_data->llg->copyTimer.start();
	thrust::copy(udata, udata + N, mag_vec_tmp.begin());
	u_data->llg->copyTimer.add();
	u_data->llg->operator()(mag_vec_tmp, ret_vec_tmp, t); // integrate
	u_data->llg->copyTimer.start();
	cudaMemcpy(dudata, thrust::raw_pointer_cast(ret_vec_tmp.data()), N * sizeof(realtype), cudaMemcpyDeviceToDevice);
	u_data->llg->copyTimer.add();
	return (0); // <--- This is required to signal success
}
