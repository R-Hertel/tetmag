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
 * TheLLG_GPU.cu
 *
 * GPU implementations of TheLLG methods and GPU_Integrator.
 * Separate translation unit required because nvcc must compile anything
 * referencing thrust::device_vector or CUDA N_Vector types.
 */

#include "LLGIntegrator.h"
#include "TheLLG.h"

#include <functional>
#include <iostream>

#include <thrust/copy.h>
#include <thrust/transform.h>
#include <thrust/functional.h>
#include <thrust/execution_policy.h>
#include <nvector/nvector_cuda.h>

typedef thrust::device_vector<double> dev_vec;


// =============================================================================
// TheLLG -- GPU LLG equation variants
// =============================================================================

dev_vec TheLLG::classicVersion_GPU(const dev_vec& mag_vec)
{
    gpucalc->setMagDev(mag_vec);
    gpucalc->computeAndAccumulateHeff(useUniaxial, useDMI);
    Eigen::MatrixXd Hcpu = assembleCpuOnlyFields();
    gpucalc->addHostContribution(Hcpu);
    return *gpucalc->ClassicLLG_dev(alpha);
}

dev_vec TheLLG::noPrecession_GPU(const dev_vec& mag_vec)
{
    gpucalc->setMagDev(mag_vec);
    gpucalc->computeAndAccumulateHeff(useUniaxial, useDMI);
    Eigen::MatrixXd Hcpu = assembleCpuOnlyFields();
    gpucalc->addHostContribution(Hcpu);
    return *gpucalc->LLG_noPrec_dev(alpha);
}

dev_vec TheLLG::sttDynamics_GPU(const dev_vec& mag_vec)
{
    dev_vec LLGpart_d = classicVersion_GPU(mag_vec);
    gpucalc->UTermSTT_GPU();
    dev_vec ret_vec_d = *gpucalc->STT_term_LLG_dev(alpha, stt.beta);
    double pulseVal = stt.sttPulse ? stt.gaussPulseValue(timeInPs) : 0.0;
    double currentVal = stt.dcAmplitude + stt.pulseAmplitude * pulseVal;
    thrust::transform(ret_vec_d.begin(), ret_vec_d.end(),
                      ret_vec_d.begin(),
                      currentVal * thrust::placeholders::_1);
    thrust::transform(thrust::device,
                      ret_vec_d.begin(), ret_vec_d.end(),
                      LLGpart_d.begin(), ret_vec_d.begin(),
                      thrust::plus<double>());
    return ret_vec_d;
}

Eigen::MatrixXd TheLLG::assembleCpuOnlyFields()
{
    if (hp.pulseIsUsed) calcPulseField();
    if (hp.sweepIsUsed) calcSweepField();
    if (hp.rfIsUsed)    calcRFField();
    return Hext + Hdem + Hpls + Hswp + Hstat + Hrf;
}

void TheLLG::selectLLGTypeGPU(int choice)
{
    enum choices { noPrec, usualLLG, STT };
    useSTT = false;
    if (choice == noPrec) {
        selectedLLGType_GPU = [this](const dev_vec& m) -> dev_vec {
            return noPrecession_GPU(m);
        };
    } else if (choice == STT) {
        useSTT = true;
        selectedLLGType_GPU = [this](const dev_vec& m) -> dev_vec {
            return sttDynamics_GPU(m);
        };
    } else {
        selectedLLGType_GPU = [this](const dev_vec& m) -> dev_vec {
            return classicVersion_GPU(m);
        };
    }
}

void TheLLG::operator()(const dev_vec& mag_d, dev_vec& dxdt_d, const double theTime)
{
    timeInPs = theTime / realTimeScale;
    dxdt_d   = selectedLLGType_GPU(mag_d);
}


void TheLLG::updateDeviceMag(const dev_vec& mag_d)
{
    gpucalc->setMagDev(mag_d);
}

// =============================================================================
// GPU_Integrator
// =============================================================================

GPU_Integrator::GPU_Integrator(TheLLG* llg, int nx)
    : llg_(llg), nx_(nx), m_gpu_(NULL), cvode_mem_(NULL), LS_(NULL), data_(), sunctx_(NULL)
{
    mag_vec_tmp_.resize(3 * nx);
    ret_vec_tmp_.resize(3 * nx);
}

GPU_Integrator::~GPU_Integrator()
{
    if (m_gpu_    != NULL) N_VDestroy_Cuda(m_gpu_);
    if (cvode_mem_!= NULL) CVodeFree(&cvode_mem_);
    if (LS_       != NULL) SUNLinSolFree(LS_);
    if (sunctx_   != NULL) SUNContext_Free(&sunctx_);
}

// GPUUserData extends UserData to carry a GPU_Integrator* so rhs_d can
// access mag_vec_tmp_ and ret_vec_tmp_ without file-scope globals.
struct GPUUserData : public UserData {
    GPU_Integrator* integrator;
};

std::shared_ptr<UserData> GPU_Integrator::allocUserData(int nx)
{
    std::shared_ptr<GPUUserData> data = std::make_shared<GPUUserData>();
    data->llg        = llg_;
    data->nx         = nx;
    data->integrator = this;
    data->ret.resize(3 * nx);
    data->mag.resize(3 * nx);
    return data;
}

int GPU_Integrator::rhs_d(realtype t, N_Vector u, N_Vector u_dot, void* user_data)
{
    (void) t;
    GPUUserData* gpu_data = static_cast<GPUUserData*>(user_data);
    GPU_Integrator* self  = gpu_data->integrator;
    sunindextype N        = 3 * gpu_data->nx;

    realtype* udata  = N_VGetDeviceArrayPointer_Cuda(u);
    realtype* dudata = N_VGetDeviceArrayPointer_Cuda(u_dot);

    thrust::copy_n(udata, N, self->mag_vec_tmp_.begin());
    gpu_data->llg->operator()(self->mag_vec_tmp_, self->ret_vec_tmp_, t);
    cudaMemcpy(dudata,
               thrust::raw_pointer_cast(self->ret_vec_tmp_.data()),
               N * sizeof(realtype),
               cudaMemcpyDeviceToDevice);
    return 0;
}

int GPU_Integrator::integrateCVODE(std::vector<double>& mag_vec,
                                    double ode_start_t,
                                    double ode_end_t,
                                    double dt)
{
    (void) dt;
    long int its_l  = 0;
    long int its_nl = 0;
    int flag;
    sunindextype N = static_cast<sunindextype>(3 * nx_);

    realtype t0     = ode_start_t;
    realtype t      = ode_start_t;
    realtype reltol = 1.0e-6;
    realtype abstol = 1.0e-6;
    realtype tout   = ode_start_t + ode_end_t;

    mag_vec_tmp_ = mag_vec;

    if (cvode_mem_ == NULL)
    {
        flag = SUNContext_Create(NULL, &sunctx_);
        if (flag != 0)
            return 1;

        m_gpu_ = N_VMake_Cuda(N, mag_vec.data(),
                               thrust::raw_pointer_cast(mag_vec_tmp_.data()),
                               sunctx_);
        if (m_gpu_ == NULL)
            return 1;

        cvode_mem_ = CVodeCreate(CV_ADAMS, sunctx_);
        if (cvode_mem_ == NULL)
            return 1;

        flag = CVodeInit(cvode_mem_, GPU_Integrator::rhs_d, t0, m_gpu_);
        if (flag != CV_SUCCESS)
            return 1;

        flag = CVodeSStolerances(cvode_mem_, reltol, abstol);
        if (flag != CV_SUCCESS)
            return 1;

        data_ = allocUserData(nx_);
        flag  = CVodeSetUserData(cvode_mem_, static_cast<void*>(data_.get()));
        if (flag != CV_SUCCESS)
            return 1;

        LS_ = SUNLinSol_SPGMR(m_gpu_, PREC_NONE, 0, sunctx_);
        if (LS_ == NULL)
            return 1;

        flag = CVodeSetLinearSolver(cvode_mem_, LS_, NULL);
        if (flag != CV_SUCCESS)
            return 1;
    }
    else
    {
        flag = CVodeReInit(cvode_mem_, t0, m_gpu_);
        if (flag != CV_SUCCESS)
            return 1;
    }

    t    = t0;
    flag = CVode(cvode_mem_, tout, m_gpu_, &t, CV_NORMAL);
    if (flag < 0) {
        std::cerr << "Warning: GPU integration failed, flag = " << flag << std::endl;
        return 1;
    }

    CVodeGetNumNonlinSolvIters(cvode_mem_, &its_nl);
    CVodeGetNumLinIters(cvode_mem_, &its_l);

    realtype* res_p = N_VGetDeviceArrayPointer_Cuda(m_gpu_);
    thrust::copy(res_p, res_p + N, mag_vec_tmp_.begin());
    thrust::copy(mag_vec_tmp_.begin(), mag_vec_tmp_.end(), mag_vec.begin());

    llg_->updateDeviceMag(mag_vec_tmp_);

    int its = static_cast<int>(its_nl + its_l);
    return its;
}
