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
 * LLGIntegrator.h
 */

#ifndef LLGINTEGRATOR_H_
#define LLGINTEGRATOR_H_

#include <vector>
#include <memory>
#include <functional>

#include <sunlinsol/sunlinsol_spgmr.h>
#include <cvode/cvode_spils.h>
#include <cvode/cvode.h>
#include <nvector/nvector_serial.h>
#include <sundials/sundials_types.h>
#include <sundials/sundials_version.h>
#include <sunlinsol/sunlinsol_dense.h>
#include <sunmatrix/sunmatrix_dense.h>
#include <nvector/nvector_openmp.h>

#ifdef USE_CUDA
#include <nvector/nvector_cuda.h>
#include <thrust/device_vector.h>
#include <thrust/device_ptr.h>
#endif

#include "typedefs.h"

typedef std::vector<double> state_type;

class TheLLG;


class UserData {
public:
    TheLLG*             llg;
    int                 nx;
    std::vector<double> ret;
    std::vector<double> mag;
};


class CPU_Integrator {
public:
    CPU_Integrator(TheLLG* llg, int nx);
    ~CPU_Integrator();

    int integrateCVODE(state_type& mag_vec,
                       double ode_start_t,
                       double ode_end_t,
                       double dt);

private:
    TheLLG*                     llg_;
    int                         nx_;

    N_Vector                    m_;
    void*                       cvode_mem_;
    SUNLinearSolver             LS_;
    std::shared_ptr<UserData>   data_;
    SUNContext                  sunctx_;

    static int rhs(realtype t, N_Vector u, N_Vector u_dot, void* user_data);

    std::shared_ptr<UserData> allocUserData(int nx);
};


#ifdef USE_CUDA
class GPU_Integrator {
public:
    GPU_Integrator(TheLLG* llg, int nx);
    ~GPU_Integrator();

    int integrateCVODE(state_type& mag_vec,
                       double ode_start_t,
                       double ode_end_t,
                       double dt);

private:
    TheLLG*                     llg_;
    int                         nx_;

    N_Vector                    m_gpu_;
    void*                       cvode_mem_;
    SUNLinearSolver             LS_;
    std::shared_ptr<UserData>   data_;
    SUNContext                  sunctx_;

    thrust::device_vector<double> mag_vec_tmp_;
    thrust::device_vector<double> rhs_mag_tmp_;
    thrust::device_vector<double> ret_vec_tmp_;

    static int rhs_d(realtype t, N_Vector u, N_Vector u_dot, void* user_data);

    std::shared_ptr<UserData> allocUserData(int nx);
};
#endif /* USE_CUDA */

#endif /* LLGINTEGRATOR_H_ */
