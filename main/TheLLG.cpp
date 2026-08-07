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
 * TheLLG.cpp
 *
 *  Created on: May 10, 2017
 *      Author: riccardo
 */

#include "TheLLG.h"
#include "LLGIntegrator.h"
#include "auxiliaries.h"
#include "PhysicalConstants.h"

#include <iostream>
#include <iomanip>
#include <Eigen/Dense>
#include <cmath>
#include <algorithm>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace Eigen;

static double odeTime;

enum Coords { x, y, z };


// =============================================================================
// Constructor / Destructor
// =============================================================================

TheLLG::~TheLLG() = default;

TheLLG::TheLLG(SimulationData& sd, const MeshData& msh, int LLGVersion)
    : nx(sd.nx),
      NodeVolume(msh.NodeVolume),
      alpha(sd.alpha),
      Js(sd.Js),
      hp(sd.Hp),
      freezeDemag(true)
{
    useGPU = sd.useGPU;
#ifdef USE_CUDA
    if (useGPU) {
        gpucalc = std::make_shared<EffFieldGPU>(nx);
    }
#endif
    efc = std::make_shared<EffFieldCalc>(sd, msh);

    odeTime      = 0.;
    refresh_time = 0.1;

    heffTimer.reset();
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
    Hstat = MatrixXd::Zero(nx, 3);
    Hrf   = MatrixXd::Zero(nx, 3);
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

    useCubic             = !(sd.Kc1.isZero() && sd.Kc2.isZero());
    useUniaxial          = !sd.Ku1.isZero();
    useSurfaceAnisotropy = !sd.Ks.isZero();

    if (useUniaxial && useGPU) {
#ifdef USE_CUDA
        gpucalc->setUniaxialAnisotropy(sd.Kuni, sd.Ku1);
#endif
    }

    realTimeScale = sd.psTimeScaleFactor();
    invJs = Js.cwiseInverse().unaryExpr([](double v) {
        return std::isfinite(v) ? v : 0.0;
    });
#ifdef USE_CUDA
    if (useGPU) {
        gpucalc->setInvJs(invJs);
    }
#endif
    JsVol = Js.cwiseProduct(NodeVolume);

    if (hp.pulseIsUsed || hp.sweepIsUsed || hp.rfIsUsed || hp.pulseHasProfile || hp.rfHasProfile)
        setHdynField();
    if (hp.staticLocalIsUsed)
        Hstat = hp.staticLocalH;

    totalVolume = NodeVolume.sum();
    selectEffectiveFields();
}


// =============================================================================
// Operator ()
// =============================================================================

void TheLLG::operator()(const state_type& state, state_type& dxdt, const double theTime)
{
    timeInPs = theTime / realTimeScale;
    odeTime  = timeInPs;
    dxdt     = selectedLLGType(state);
}


// =============================================================================
// LLG type selection
// =============================================================================

void TheLLG::selectLLGType(int choice)
{
    enum choices { noPrec, usualLLG, STT };
    useSTT = false;

    if (useGPU) {
#ifdef USE_CUDA
        selectLLGTypeGPU(choice);
#endif
    }

    if (choice == noPrec) {
        selectedLLGType = [this](const state_type& m) -> state_type {
            return noPrecession(m);
        };
    } else if (choice == STT) {
        useSTT = true;
        selectedLLGType = [this](const state_type& m) -> state_type {
            return sttDynamics(m);
        };
    } else {
        selectedLLGType = [this](const state_type& m) -> state_type {
            return classicVersion(m);
        };
    }
}


// =============================================================================
// LLG equation variants
// =============================================================================

state_type TheLLG::classicVersion(const state_type& mag_vec)
{
    Map<const MatrixXd_CM> Mag(mag_vec.data(), nx, 3);
    evaluateAllEffectiveFields(Mag);
    Heff = totalEffectiveField();
    const MatrixXd MxH = cross(Mag, Heff);
    Map<MatrixXd_CM>(ret_vec.data(), ret.rows(), 3) = -MxH - alpha * cross(Mag, MxH);
    return ret_vec;
}

state_type TheLLG::noPrecession(const state_type& mag_vec)
{
    Map<const MatrixXd_CM> Mag(mag_vec.data(), nx, 3);
    evaluateAllEffectiveFields(Mag);
    Heff = totalEffectiveField();
    Map<MatrixXd_CM>(ret_vec.data(), ret.rows(), 3) = -alpha * cross(Mag, cross(Mag, Heff));
    return ret_vec;
}

state_type TheLLG::sttDynamics(const state_type& mag_vec)
{
    const std::vector<double> LLGpart = classicVersion(mag_vec);
    Map<const MatrixXd_CM> Mag(mag_vec.data(), nx, 3);
    calcUtermSTT(Mag);
    double pulseVal = stt.sttPulse ? stt.gaussPulseValue(timeInPs) : 0.0;
    double currentVal = stt.dcAmplitude + stt.pulseAmplitude * pulseVal;
    const MatrixXd MxU = cross(Mag, stt.Ustt);
    ret = -(stt.beta - alpha) * MxU - (1. + alpha * stt.beta) * cross(Mag, MxU);
    ret *= currentVal;
    Map<MatrixXd_CM>(ret_vec.data(), ret.rows(), 3) = ret;
    std::transform(ret_vec.begin(), ret_vec.end(),
                   LLGpart.begin(), ret_vec.begin(), std::plus<double>());
    return ret_vec;
}

void TheLLG::calcUtermSTT(MRef& Mag)
{
    for (int i = 0; i < 3; ++i) {
        stt.Ustt.col(i) = stt.eta_jx.cwiseProduct(stt.gradX * Mag.col(i))
                        + stt.eta_jy.cwiseProduct(stt.gradY * Mag.col(i))
                        + stt.eta_jz.cwiseProduct(stt.gradZ * Mag.col(i));
    }
}


// =============================================================================
// Effective field assembly
// =============================================================================

void TheLLG::selectEffectiveFields()
{
    calcEffectiveField.clear();
    calcEffectiveField.emplace_back([this](MRef& m) { calcExchangeField(m); });
    if (useUniaxial)
        calcEffectiveField.emplace_back([this](MRef& m) { calcUniaxialAnisotropyField(m); });
    if (useCubic)
        calcEffectiveField.emplace_back([this](MRef& m) { calcCubicAnisotropyField(m); });
    if (useDMI)
        calcEffectiveField.emplace_back([this](MRef& m) { calcDMIField(m); });
    if (useSurfaceAnisotropy)
        calcEffectiveField.emplace_back([this](MRef& m) { calcSurfaceAnisotropyField(m); });
    if (hp.pulseIsUsed)
        calcEffectiveField.emplace_back([this](MRef&) { calcPulseField(); });
    if (hp.sweepIsUsed)
        calcEffectiveField.emplace_back([this](MRef&) { calcSweepField(); });
    if (hp.rfIsUsed)
        calcEffectiveField.emplace_back([this](MRef&) { calcRFField(); });
}

void TheLLG::evaluateAllEffectiveFields(MRef& Mag)
{
    int tasks = calcEffectiveField.size();
    for (int i = 0; i < tasks; ++i) {
        calcEffectiveField[i](Mag);
    }
}

MatrixXd TheLLG::totalEffectiveField()
{
    MatrixXd Htot(nx, 3);
#ifdef _OPENMP
    int i, j;
#pragma omp parallel for private(i) collapse(2)
    for (i = 0; i < nx; ++i) {
        for (j = 0; j < 3; ++j) {
            Htot(i, j) = (Hexc(i,j) + Hku(i,j) + Hcub(i,j) + Hsurf(i,j) + Hdmi(i,j)) * invJs(i)
                       + Hext(i,j) + Hdem(i,j) + Hpls(i,j) + Hswp(i,j) + Hstat(i,j) + Hrf(i,j);
        }
    }
#else
    MatrixXd Heff_tot = (Hexc + Hku + Hcub + Hsurf + Hdmi).array().colwise() * invJs.array();
    Htot = Heff_tot + Hext + Hdem + Hpls + Hswp + Hstat + Hrf;
#endif
    return Htot;
}

void TheLLG::calcExchangeField(MRef& Mag)
{
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

void TheLLG::calcDemagField(MRef& Mag)
{
    if (odeTime >= next_refresh) {
        next_refresh += refresh_time;
        Hdem = demag->calcField(Mag);
    }
}

void TheLLG::calcUniaxialAnisotropyField(MRef& Mag)
{
    if (useGPU) {
#ifdef USE_CUDA
        Hku = gpucalc->UniaxialAnisotropyField();
#endif
    } else {
        Hku = efc->uniaxialAnisotropyField(Mag);
    }
}

void TheLLG::calcCubicAnisotropyField(MRef& Mag)
{
    if (useCubic) {
        Hcub = efc->cubicAnisotropyField(Mag);
    }
}

void TheLLG::calcSurfaceAnisotropyField(MRef& Mag)
{
    Hsurf = efc->surfaceAnisotropyField(Mag);
}

void TheLLG::calcDMIField(MRef& Mag)
{
    if (useGPU) {
#ifdef USE_CUDA
        Hdmi = gpucalc->DMIField();
#endif
    } else {
        Hdmi = efc->dmiField(Mag);
    }
}

void TheLLG::calcPulseField()
{
    Hpls = hp.pulseH * hp.gaussPulseValue(timeInPs);
}

void TheLLG::calcSweepField()
{
    Hswp = hp.sweepH * hp.sweepFieldValue(timeInPs);
}

void TheLLG::calcRFField()
{
    Hrf = hp.rfH * hp.rfFieldValue(timeInPs);
}


// =============================================================================
// Energies
// =============================================================================

double TheLLG::getExchangeEnergy(MRef& Mag)
{
    return efc->calcExchangeEnergy(Mag);
}

double TheLLG::getDirectExch(MRef& Mag)
{
    return efc->exchEnergy_direct(Mag);
}

double TheLLG::getZeemanEnergy(MRef& Mag)
{
    double zeeEnergy = -(JsVol.transpose() * (Hext + Hpls + Hswp + Hstat + Hrf).cwiseProduct(Mag)).sum();
    return zeeEnergy;
}

double TheLLG::getDemagEnergy(MRef& Mag)
{
    if (freezeDemag)
        calcDemagField(Mag);
    return demag->getDemagEnergy(Mag);
}

double TheLLG::getUniaxialAnisotropyEnergy(MRef& Mag)
{
    if (!useUniaxial)
        return 0;
    return efc->uniaxialEnergy(Mag);
}

double TheLLG::getCubicAnisotropyEnergy(MRef& Mag)
{
    if (!useCubic)
        return 0;
    return efc->cubicAnisotropyEnergy(Mag);
}

double TheLLG::getDMIEnergy(MRef& Mag)
{
    if (!useDMI)
        return 0;
#ifdef USE_CUDA
    calcDMIField(Mag);
    efc->setDMIField(Hdmi);
#endif
    return efc->dmiEnergy(Mag);
}

double TheLLG::getSurfaceAnisotropyEnergy(MRef& Mag)
{
    if (!useSurfaceAnisotropy)
        return 0.;
    return efc->surfaceAnisotropyEnergy(Mag);
}

double TheLLG::getDirectDMI(MRef& Mag)
{
    return efc->dmiEnergy_direct(Mag);
}


// =============================================================================
// Getters
// =============================================================================

Eigen::Vector3d TheLLG::getMeanH()
{
    return NodeVolume.transpose() * (Hext + Hpls + Hswp + Hstat + Hrf) / totalVolume;
}

double TheLLG::getMaxTorque(MRef& Mag)
{
    evaluateAllEffectiveFields(Mag);
    Heff = totalEffectiveField();
    double m_torque = 0;
    if (useGPU) {
#ifdef USE_CUDA
        m_torque = gpucalc->MaxTorque(Heff);
#endif
    } else {
        m_torque = cross(Mag, Heff).rowwise().norm().maxCoeff();
    }
    return m_torque * PhysicalConstants::mu0;
}

void TheLLG::outputTimer()
{
    std::cout << "time for copying data [s]:\t"
              << std::setprecision(2) << copyTimer.durationInMus() / 1.e6 << std::endl;
#ifdef USE_CUDA
    gpucalc->displayTimer();
#endif
    std::cout << "time for effective field: "
              << std::setprecision(2) << heffTimer.durationInMus() / 1.e6 << std::endl;
    std::cout << std::setprecision(-1);
}


// =============================================================================
// Setters
// =============================================================================

void TheLLG::setTime(double time)
{
    timeInPs = time;
}

void TheLLG::setDemagData(std::shared_ptr<DemagField> demag_)
{
    demag = demag_;
    freezeDemag = false;
    calcEffectiveField.emplace_back([this](MRef& m) { calcDemagField(m); });
}

void TheLLG::setSTTData(const STT& stt_)
{
    stt = stt_;
    stt.setEta(invJs);
    if (useGPU) {
#ifdef USE_CUDA
        std::cout << "setting STT data on device" << std::endl;
        gpucalc->setSTTDataOnDevice(stt.gradX, stt.gradY, stt.gradZ,
                                    stt.eta_jx, stt.eta_jy, stt.eta_jz);
#endif
    }
}

namespace {
Eigen::MatrixXd fieldFromAngles(int nx, double Habs, double theta_h, double phi_h)
{
    MatrixXd Hfield = MatrixXd::Ones(nx, 3);
    Hfield.col(x) *= std::sin(theta_h) * std::cos(phi_h);
    Hfield.col(y) *= std::sin(theta_h) * std::sin(phi_h);
    Hfield.col(z) *= std::cos(theta_h);
    Hfield *= (Habs / PhysicalConstants::mu0);
    return Hfield;
}
} // anonymous namespace

void TheLLG::setZeemanField(int nx, double H_stat, double theta_H_deg, double phi_H_deg,
                             double H_hys,  double theta_Hys_deg, double phi_Hys_deg)
{
    double deg_to_rad = PhysicalConstants::pi / 180.;
    MatrixXd Hhys  = fieldFromAngles(nx, H_hys,  theta_Hys_deg * deg_to_rad, phi_Hys_deg * deg_to_rad);
    MatrixXd Hstat = fieldFromAngles(nx, H_stat, theta_H_deg   * deg_to_rad, phi_H_deg   * deg_to_rad);
    Hext = Hhys + Hstat;
}

void TheLLG::setHdynField()
{
    hp.pulseH = MatrixXd::Zero(nx, 3);
    hp.sweepH = MatrixXd::Zero(nx, 3);
    hp.rfH    = MatrixXd::Zero(nx, 3);
    if (hp.pulseIsUsed) {
        if (hp.pulseHasProfile)
            hp.pulseH = hp.fieldProfile / PhysicalConstants::mu0;
        else
            hp.pulseH = fieldFromAngles(nx, 1., hp.pulseTheta, hp.pulsePhi);
    }
    if (hp.sweepIsUsed)
        hp.sweepH = fieldFromAngles(nx, 1., hp.sweepTheta, hp.sweepPhi);
    if (hp.rfIsUsed) {
        if (hp.rfHasProfile)
            hp.rfH = hp.fieldProfile / PhysicalConstants::mu0;
        else
            hp.rfH = fieldFromAngles(nx, 1., hp.rfTheta, hp.rfPhi);
    }
}

void TheLLG::setMag(MRef& mag_)
{
    (void) mag_;
    if (useGPU) {
#ifdef USE_CUDA
        gpucalc->setMagDev(mag_);
#endif
    }
}

void TheLLG::setHdem(MRef& Hdem_)
{
    Hdem = Hdem_;
}

// =============================================================================
// Integrator management
// =============================================================================

void TheLLG::initIntegrator()
{
#ifdef USE_CUDA
    if (useGPU) {
        gpu_integrator_ = std::make_unique<GPU_Integrator>(this, nx);
        integrator_ = [this](std::vector<double>& mag_vec,
                              double ode_start_t, double ode_end_t, double dt) -> int {
            return gpu_integrator_->integrateCVODE(mag_vec, ode_start_t, ode_end_t, dt);
        };
        return;
    }
#endif
    cpu_integrator_ = std::make_unique<CPU_Integrator>(this, nx);
    integrator_ = [this](std::vector<double>& mag_vec,
                          double ode_start_t, double ode_end_t, double dt) -> int {
        return cpu_integrator_->integrateCVODE(mag_vec, ode_start_t, ode_end_t, dt);
    };
}

int TheLLG::integrateSUNDIALS(std::vector<double>& mag_vec,
                                double ode_start_t, double ode_end_t, double dt)
{
    copyTimer.start();
    int its = integrator_(mag_vec, ode_start_t, ode_end_t, dt);
    copyTimer.add();
    return its;
}
