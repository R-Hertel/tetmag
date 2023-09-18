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
#include <functional> // for std::function, ref, reference_wrapper
#include "SimulationData.h"
#include "MeshData.h"
#include "typedefs.h"
#ifdef USE_CUDA 
	#include <memory>
	#include "EffFieldGPU.h"
#endif
#include "DemagField.h"

#include <sunlinsol/sunlinsol_spgmr.h>  //access to SPGMR SUNLinearSolver
#include <cvode/cvode_spils.h> // access to CVSpils interface
#include <cvode/cvode.h> // prototypes for CVODE fcts., consts.
#include <nvector/nvector_serial.h>  // access to serial N_Vector
#include <sundials/sundials_version.h>  // defs. of realtype, sunindextype
#include <sundials/sundials_types.h>
#include <sunlinsol/sunlinsol_dense.h> /* access to dense SUNLinearSolver          */
#include <sunmatrix/sunmatrix_dense.h> /* access to dense SUNMatrix                */
#include <nvector/nvector_openmp.h>    /* serial N_Vector types, fcts., macros */

#include "Timer.h"
#include <memory>



struct UserData;

typedef std::vector<double> state_type;
#include <EffFieldCalc.h>

class TheLLG : public std::enable_shared_from_this<TheLLG> {
private:
	Timer copyTimer, heffTimer;
	int nx;
	double next_refresh;
	double refresh_time;
	Eigen::VectorXd NodeVolume;
	double alpha;
	double Habs;
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
	Eigen::MatrixXd Hloc;
	double theta_H, phi_H;
	Eigen::VectorXd Js;
	Eigen::VectorXd JsVol;
	Eigen::VectorXd invNodeVol;
	Eigen::VectorXd invJs;
	Eigen::MatrixXd ret;
	std::vector<double> ret_vec;
	Hdynamic hp;
	Hlocal hl;
	std::shared_ptr<EffFieldCalc> efc;
	state_type classicVersion(const state_type&);
	state_type noPrecession(const state_type&);
	state_type sttDynamics(const state_type&);
	void calcUtermSTT(MRef&);
#ifdef USE_CUDA
	std::shared_ptr<EffFieldGPU> gpucalc;
	typedef thrust::device_vector<double> dev_vec;
	dev_vec classicVersion_GPU(const dev_vec &);
	dev_vec noPrecession_GPU(const dev_vec &);
	dev_vec sttDynamics_GPU(const dev_vec &);
	std::function <dev_vec(const dev_vec&)> selectedLLGType_GPU;
	Eigen::MatrixXd effFieldsForGPU(const dev_vec&);
	static int rhs_d(realtype t, N_Vector u, N_Vector u_dot, void *user_data);
#endif
	std::function <state_type(const state_type&)> selectedLLGType;
	std::vector <std::function <void(MRef&)> > calcEffectiveField;
	void selectLLGType(int);
	void calcExchangeField(MRef&);
	bool useUniaxial;
	void calcUniaxialAnisotropyField(MRef&);
	bool useCubic;
	void calcCubicAnisotropyField(MRef&);
	bool useSurfaceAnisotropy;
	void calcSurfaceAnisotropyField(MRef&);
	bool useDMI;
	bool useGPU;
	bool noDemag;
	void calcDMIField(MRef&);
	void calcPulseField(MRef&);
	void calcSweepField(MRef&);
	void calcLocalField(MRef&);
	double totalVolume;
	STT stt; // struct defined in ProgramSpecs.h
	std::shared_ptr<DemagField> demag;
	bool freezeDemag;
	void calcDemagField(MRef&);
	void evaluateAllEffectiveFields(MRef&);
	Eigen::MatrixXd totalEffectiveField();
	double realTimeScale;
	double timeInPs;
	void setHdynField();
	bool useSTT;
	void selectEffectiveFields();
	void selectLLGTypeGPU(int);
	static int rhs(realtype t, N_Vector u, N_Vector u_dot, void *user_data);
//	std::shared_ptr<TheLLG> getptr();
//	UserData* alloc_user_data(const int );
	std::shared_ptr<UserData> alloc_user_data(const int );
public:
//	static TheLLG* llg_p; // unused, attempt to transfer 'this' to static member function in CVODE
	TheLLG(SimulationData& , const MeshData&, int );
	void operator()(const state_type&, state_type&, const double /*time*/);
#ifdef USE_CUDA 
    void operator()(const dev_vec&, dev_vec&, const double /*time*/);
#endif
    int IntegrateOnGPU( std::vector<double>& , double , double , double );
    int gpuODE( std::vector<double>& , double , double , double );
    int IntegrateODEINT( std::vector<double>& , double , double , double );
    int integrateSUNDIALS( std::vector<double>& , double , double , double );
	double getDemagEnergy(MRef&);
	double getDMIEnergy(MRef&);
	double getSurfaceAnisotropyEnergy(MRef&);
	double getCubicAnisotropyEnergy(MRef&);
	double getUniaxialAnisotropyEnergy(MRef&);
	double getExchangeEnergy(MRef&);
	double getZeemanEnergy(MRef&);
	double getMaxTorque(MRef&);
	void setDemagData( const DemagField& );
	void setSTTData( const STT& );
	void setZeemanField(int, double, double, double, double, double, double);
	void setHdem(MRef&);
	void setTime(double);
	void setMag(MRef&);
	Eigen::Vector3d getMeanH();
	void outputTimer();
	double getDirectExch(MRef&);
	double getDirectDMI(MRef&);
};


struct UserData{
  std::shared_ptr<TheLLG> llg;
  int nx;
  std::vector<double> ret;
};


#endif /* THELLG_H_ */
