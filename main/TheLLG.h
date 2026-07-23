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
 * TheLLG.h
 *
 *  Created on: May 10, 2017
 *      Author: riccardo
 */

#ifndef THELLG_H_
#define THELLG_H_

#include <Eigen/SparseCore>
#include <Eigen/Dense>
#include <vector>
#include <functional>
#include <memory>

#include "SimulationData.h"
#include "MeshData.h"
#include "typedefs.h"
#include "DemagField.h"
#include "EffFieldCalc.h"
#include "Timer.h"

#ifdef USE_CUDA
#include "EffFieldGPU.h"
#endif

class CPU_Integrator;
#ifdef USE_CUDA
class GPU_Integrator;
#endif

typedef std::vector<double> state_type;


class TheLLG : public std::enable_shared_from_this<TheLLG> {
public:
    TheLLG(SimulationData&, const MeshData&, int);
    ~TheLLG();

    void operator()(const state_type&, state_type&, const double /*time*/);
#ifdef USE_CUDA
    typedef thrust::device_vector<double> dev_vec;
    void operator()(const dev_vec&, dev_vec&, const double /*time*/);
#endif

    void initIntegrator();
    int  integrateSUNDIALS(std::vector<double>&, double, double, double);

    double getDemagEnergy(MRef&);
    double getDMIEnergy(MRef&);
    double getSurfaceAnisotropyEnergy(MRef&);
    double getCubicAnisotropyEnergy(MRef&);
    double getUniaxialAnisotropyEnergy(MRef&);
    double getExchangeEnergy(MRef&);
    double getZeemanEnergy(MRef&);
    double getDirectExch(MRef&);
    double getDirectDMI(MRef&);
    double getMaxTorque(MRef&);

    void setDemagData(const DemagField&);
    void setSTTData(const STT&);
    void setZeemanField(int, double, double, double, double, double, double);
    void setHdem(MRef&);
    void setTime(double);
    void setMag(MRef&);
#ifdef USE_CUDA
    void updateDeviceMag(const thrust::device_vector<double>& mag_d);
#endif

    Eigen::Vector3d getMeanH();
    void outputTimer();

private:
    Timer copyTimer, heffTimer;

    int                 nx;
    Eigen::VectorXd     NodeVolume;
    double              alpha;
    Eigen::VectorXd     Js;
    Eigen::VectorXd     JsVol;
    Eigen::VectorXd     invNodeVol;
    Eigen::VectorXd     invJs;
    double              totalVolume;
    double              realTimeScale;
    double              timeInPs;

    Eigen::MatrixXd Hexc;
    Eigen::MatrixXd Hdem;
    Eigen::MatrixXd Hku;
    Eigen::MatrixXd Hext;
    Eigen::MatrixXd Hcub;
    Eigen::MatrixXd Hsurf;
    Eigen::MatrixXd Heff;
    Eigen::MatrixXd Hdmi;
    Eigen::MatrixXd Hpls;
    Eigen::MatrixXd Hswp;
    Eigen::MatrixXd Hstat;
    Eigen::MatrixXd Hrf;

    Hdynamic hp;

    std::shared_ptr<EffFieldCalc>  efc;
    std::shared_ptr<DemagField>    demag;
    bool                           freezeDemag;
    double                         next_refresh;
    double                         refresh_time;

    bool useUniaxial;
    bool useCubic;
    bool useSurfaceAnisotropy;
    bool useDMI;
    bool useSTT;
    bool useGPU;

    STT stt;

    Eigen::MatrixXd     ret;
    std::vector<double> ret_vec;

    state_type classicVersion(const state_type&);
    state_type noPrecession(const state_type&);
    state_type sttDynamics(const state_type&);
    void       calcUtermSTT(MRef&);

    std::function<state_type(const state_type&)> selectedLLGType;

    void selectLLGType(int);

    std::vector<std::function<void(MRef&)>> calcEffectiveField;

    void selectEffectiveFields();
    void evaluateAllEffectiveFields(MRef&);
    Eigen::MatrixXd totalEffectiveField();

    void calcExchangeField(MRef&);
    void calcDemagField(MRef&);
    void calcUniaxialAnisotropyField(MRef&);
    void calcCubicAnisotropyField(MRef&);
    void calcSurfaceAnisotropyField(MRef&);
    void calcDMIField(MRef&);
    void calcPulseField();
    void calcSweepField();
    void calcRFField();

    void setHdynField();

    typedef std::function<int(std::vector<double>&, double, double, double)> IntegratorFunction;
    IntegratorFunction integrator_;

    std::unique_ptr<CPU_Integrator> cpu_integrator_;
#ifdef USE_CUDA
    std::unique_ptr<GPU_Integrator> gpu_integrator_;
#endif

#ifdef USE_CUDA
    std::shared_ptr<EffFieldGPU>    gpucalc;

    dev_vec classicVersion_GPU(const dev_vec&);
    dev_vec noPrecession_GPU(const dev_vec&);
    dev_vec sttDynamics_GPU(const dev_vec&);

    std::function<dev_vec(const dev_vec&)> selectedLLGType_GPU;

    Eigen::MatrixXd assembleCpuOnlyFields(const dev_vec&);
    void selectLLGTypeGPU(int);
#endif
};

#endif /* THELLG_H_ */
