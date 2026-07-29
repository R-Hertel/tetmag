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
 * LLGIntegrator.cpp
 */

#include "LLGIntegrator.h"
#include "TheLLG.h"

#include <algorithm>
#include <iostream>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif


CPU_Integrator::CPU_Integrator(TheLLG* llg, int nx)
    : llg_(llg), nx_(nx), m_(NULL), cvode_mem_(NULL), LS_(NULL), data_(), sunctx_(NULL)
{
}

CPU_Integrator::~CPU_Integrator()
{
#ifdef _OPENMP
    if (m_ != NULL)
        N_VDestroy_OpenMP(m_);
#else
    if (m_ != NULL)
        N_VDestroy(m_);
#endif

    if (cvode_mem_ != NULL)
        CVodeFree(&cvode_mem_);
    if (LS_ != NULL)
        SUNLinSolFree(LS_);
    if (sunctx_ != NULL)
        SUNContext_Free(&sunctx_);
}


std::shared_ptr<UserData> CPU_Integrator::allocUserData(int nx)
{
    std::shared_ptr<UserData> data = std::make_shared<UserData>();
    data->llg = llg_;
    data->nx  = nx;
    data->ret.resize(3 * nx);
    data->mag.resize(3 * nx);
    return data;
}


int CPU_Integrator::rhs(realtype t, N_Vector u, N_Vector u_dot, void* user_data)
{
    (void) t;
    UserData* u_data = static_cast<UserData*>(user_data);
    int n = 3 * u_data->nx;

    realtype* udata;
#ifdef _OPENMP
    udata = NV_DATA_OMP(u);
#else
    udata = N_VGetArrayPointer(u);
#endif

    std::vector<double>& mag = u_data->mag;
    std::copy(udata, udata + n, mag.begin());
    u_data->llg->operator()(mag, u_data->ret, t);

    realtype* dudata;
#ifdef _OPENMP
    dudata = NV_DATA_OMP(u_dot);
#else
    dudata = N_VGetArrayPointer(u_dot);
#endif

    std::copy(u_data->ret.begin(), u_data->ret.begin() + n, dudata);
    return 0;
}


int CPU_Integrator::integrateCVODE(std::vector<double>& mag_vec,
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

    long int maxNumSteps = 5000;

    if (cvode_mem_ == NULL)
    {
        flag = SUNContext_Create(NULL, &sunctx_);
        if (flag != 0)
            return 1;

#ifdef _OPENMP
        int num_threads = omp_get_max_threads();
        m_ = N_VMake_OpenMP(N, mag_vec.data(), num_threads, sunctx_);
#else
        m_ = N_VMake_Serial(N, mag_vec.data(), sunctx_);
#endif
        if (m_ == NULL)
            return 1;

        cvode_mem_ = CVodeCreate(CV_ADAMS, sunctx_);
        if (cvode_mem_ == NULL)
            return 1;

        flag = CVodeInit(cvode_mem_, CPU_Integrator::rhs, t0, m_);
        if (flag != CV_SUCCESS)
            return 1;

        flag = CVodeSStolerances(cvode_mem_, reltol, abstol);
        if (flag != CV_SUCCESS)
            return 1;

        flag = CVodeSetMaxNumSteps(cvode_mem_, maxNumSteps);
        if (flag != CV_SUCCESS)
            return 1;

        data_ = allocUserData(nx_);
        flag  = CVodeSetUserData(cvode_mem_, static_cast<void*>(data_.get()));
        if (flag != CV_SUCCESS)
            return 1;

        LS_ = SUNLinSol_SPGMR(m_, PREC_NONE, 0, sunctx_);
        if (LS_ == NULL)
            return 1;

        flag = CVodeSetLinearSolver(cvode_mem_, LS_, NULL);
        if (flag != CV_SUCCESS)
            return 1;
    }
    else
    {
        flag = CVodeReInit(cvode_mem_, t0, m_);
        if (flag != CV_SUCCESS)
            return 1;
    }

    flag = CVode(cvode_mem_, tout, m_, &t, CV_NORMAL);
    if (flag < 0) {
        std::cerr << "Warning: CPU integration reached t = " << t
                  << " instead of " << tout << ", flag = " << flag << std::endl;
        return 1;
    }

    CVodeGetNumNonlinSolvIters(cvode_mem_, &its_nl);
    CVodeGetNumLinIters(cvode_mem_, &its_l);

    int its = static_cast<int>(its_nl + its_l);
    return its;
}

// GPU_Integrator is implemented in gpu/TheLLG_GPU.cu (requires nvcc).
