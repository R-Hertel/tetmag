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
 * TheSimulation.cpp
 *
 *  Created on: Nov 30, 2016
 *      Author: riccardo
 */
#include "TheSimulation.h"
#include "InitMag.h"
#include "SimulationData.h"
#include "MeshData.h"
#include "DemagField.h"
#include "TheLLG.h"
#include "MeshWriter.h"
#include "PhysicalConstants.h"
#include "ProcessInformation.h"
#include <iostream>
#include <fstream>
#include <cmath>
#include "ProgramSpecs.h"
#include <functional>
#include "auxiliaries.h"
#include "Timer.h"
#include "Hysteresis.h"
//#include <cuda_profiler_api.h>

enum coords {x, y, z};
using namespace Eigen;

std::streamsize defaultPrecision = std::cout.precision();
std::streamsize ss = std::cout.precision();

static double realTimeScale;

TheSimulation::TheSimulation(SimulationData& sd_, MeshData& msh_, ProgramSpecs& prog_)
: nx(msh_.xyz.rows()), sd(sd_), msh(msh_), prog(prog_), totalVolume(msh.NodeVolume.sum()) {
#ifdef USE_CUDA 
	if (sd.useGPU) {
		onGPU = std::make_shared<EffFieldGPU>(nx);
	}
#endif
	DemagEnergy    = 0.;
	ExchangeEnergy = 0.;
	UniaxialEnergy = 0.;
	ZeemanEnergy   = 0.;
	CubicEnergy	   = 0.;
	SurfaceEnergy  = 0.;
	DMIEnergy 	   = 0.;
	totalEnergy    = 0.;
/*		
	directDMI	   = 0.;
	directExch	   = 0.;
*/	
}

void writeProcessInformationFile(ProcessInformation& process, std::string name) {
	std::ofstream pid(name + ".pid");
		pid << "process ID: " << process.id << "\nhost : " << process.host << "\nstarted " << process.startTime;
		pid.close();
}

void writeLogFileHeader(std::ofstream& logstream, ProcessInformation& process) {
	logstream << "# Log file for simulation started " << process.startTime;
	logstream << "# on host " << process.host << " by user " << process.user << "\n";
	logstream << "# data in the columns is: \n#(1) time in ps, (2) total energy , "
			"(3) demagnetizing energy, (4) exchange energy\n#(5) uniaxial anisotropy energy, "
			"(6) Zeeman energy, (7) cubic anisotropy energy (8) surface anisotropy energy, "
			"\n#(9) bulk DMI energy (10) maximum torque, (11-13) Mx, My, Mz, (14-16) Hx, Hy, Hz. All energies in J/m3, all fields in T"
			<< std::endl;
}

void writeLogFileFooter(std::ofstream& logstream) {
	ProcessInformation ended;
	logstream << "# Simulation ended " << ended.startTime;
}

Vector3d TheSimulation::getMeanM() {
	return msh.NodeVolume.transpose() * mag / totalVolume;
}

void myobs(std::vector<double>, double t) {
	std::cout << "t = " << t / realTimeScale << " ps " << std::endl;
}

bool timeReached(double tnow, double tschedule) {
	return (tnow > tschedule || areEqual(tnow, tschedule, 1e5));
}

void TheSimulation::start() {
//	bool integrateOnGPU = prog.integrateOnGPU;
	bool integrateOnGPU = false;
	(void) integrateOnGPU;
	if (sd.useGPU) {
		integrateOnGPU = prog.integrateOnGPU;
//		integrateOnGPU = true;
#ifdef USE_CUDA 
//	if(prog.solverType == "cuda") {
//			std::cout << "This GPU version of tetmag is powered by AMGCL." << std::endl;
//	}

#else
		std::cout << "This version of tetmag does not support GPU acceleration." << std::endl;
		exit(0);
#endif
	}
	std::cout << (sd.useGPU ? "" : "not ") << "using GPU acceleration."	<< std::endl;
	ProcessInformation process;
	writeProcessInformationFile(process, msh.name);
	realTimeScale = sd.psTimeScaleFactor();
	DemagField demag(msh, sd.Js, prog.useH2);
	std::cout << "Preparing solvers... " << std::flush;
	demag.initializeSolvers(prog.solverType, prog.preconditionerType, prog.cgTol);
	std::cout << "done." << std::endl;
//
//	bool calcDemag = (!prog.noDemag && prog.freezeDemag);
//
	bool calcDemag = !prog.noDemag;
	MeshWriter write(msh.fel, msh.xyz, msh.name, msh.NodeMaterial);

// check hysteresis (needed before calling LLG constructor)
	Hysteresis hys; // Need object even if not used.
	if (prog.useHyst) {
		hys.init(prog.firstField, prog.lastField, prog.deltaH, prog.hystTheta,	prog.hystPhi);
		writeLogFileHeader(hys.hystLogFile, process);
		std::cout << "Calculating hysteresis branch." << std::endl;
		sd.Hys = true;
		sd.H_Hys = prog.firstField / 1000.; // convert mT to T
		sd.theta_Hys = prog.hystTheta;
		sd.phi_Hys = prog.hystPhi;
		if (prog.sweepUsed)
			hys.errorMessage("sweep field");
		if (prog.pulsedSTT)
			hys.errorMessage("current type = pulse");
		if (prog.fieldPulse)
			hys.errorMessage("field pulse");
		if (prog.hl.type == "pulse")
			hys.errorMessage("local field = pulse");
		if (prog.hl.type == "sine")
			hys.errorMessage("local field = sine");
	}
	//

	write.selectWriter(prog.writerType);
	std::vector<double> mag_vec(3 * nx);
	enum choice {
		noPrecession, usualLLG, stt
	};
	int LLGType;
	STT stt_data(msh.gradX, msh.gradY, msh.gradZ);
	stt_data.setValues(sd.nx, sd.scale, prog);
	LLGType = usualLLG; // this is the default
	if (prog.noPrecession)
		LLGType = noPrecession;
	if (prog.useSTT)
		LLGType = stt;
	TheLLG LLG(sd, msh, LLGType);
	if (prog.useSTT) {
		LLG.setSTTData(stt_data);
	}
	if (!prog.freezeDemag) {
		LLG.setDemagData(demag);
	}

#ifdef USE_CUDA 
	if (sd.useGPU && integrateOnGPU) { gpuWrap.init(sd.nx); }
#endif
	LLG.setMag(mag);
	if (calcDemag) {
		LLG.setHdem(demag.calcField(mag));
	} else {
		LLG.setHdem(MatrixXd::Zero(nx, 3));
	}

	LLG.setZeemanField(nx, sd.Hext, sd.theta_H, sd.phi_H, sd.H_Hys,	sd.theta_Hys, sd.phi_Hys);

	double elapsed_t = prog.initialTime * realTimeScale;
	double elapsed_ps = elapsed_t / realTimeScale;
	double nextGraphicsOutput = elapsed_ps;
	double nextConsoleOutput = elapsed_ps;
	double nextLogOutput = elapsed_ps;
	double ode_start_t = elapsed_ps;
	const double ode_end_t = prog.timeStep * realTimeScale;
	double dt = std::min(ode_end_t, realTimeScale); // internal time step, unimportant (adaptive time steppers handle this automatically)

	std::ofstream logstream;
	if (prog.resuming) {
		logstream.open(msh.name + ".log", std::ios_base::app);
		nextGraphicsOutput += prog.cfgStride;
		nextLogOutput += prog.logStride;
	} else {
		logstream.open(msh.name + ".log");
		writeLogFileHeader(logstream, process);
	}
	logstream << std::setprecision(12);

	int fileNumber = prog.firstFileNumber;
	size_t total_its = 0;
	normalizeMag(mag);
	double maximumTorque = LLG.getMaxTorque(mag);

	if (sd.Hp.sweepIsUsed && prog.duration < sd.Hp.sweepDuration)
		prog.duration = sd.Hp.sweepDuration;
//	Timer simTimer;
	simTimer.start();
//	Timer odeTimer, pdeTimer, outputTimer, totalTimer;
	odeTimer.reset(); pdeTimer.reset(); outputTimer.reset(); totalTimer.reset();
	do {

///////////////////////////// Main Loop //////////////////////////////////////////////////////////////////
		while (!timeReached(elapsed_ps, (prog.duration + prog.initialTime))
				&& (maximumTorque > prog.maxTorque || prog.useSTT)) { // STT torque not included in maxtorque
			totalTimer.start();
			pdeTimer.start();
//		cudaProfilerStart();
//		if (prog.freezeDemag) {  // ToDo: ensure first calculation in the case of !frezee Demag
			if (calcDemag)
				LLG.setHdem(demag.calcField(mag));
//		}
		    pdeTimer.add();
			outputTimer.start();
			LLG.setTime(elapsed_ps);

/// output during simulation
			if (timeReached(elapsed_ps, nextGraphicsOutput)) {
				write.graphicsOutput(fileNumber, mag, elapsed_t / realTimeScale);
				fileNumber++;
				nextGraphicsOutput += prog.cfgStride;
			}
			if (timeReached(elapsed_ps, nextConsoleOutput)) {
				calculateEnergies(std::ref(LLG), std::ref(demag), calcDemag);
				displayEnergies(elapsed_ps);
				nextConsoleOutput += prog.consoleStride;
			}
			if (timeReached(elapsed_ps, nextLogOutput)) {
				calculateEnergies(std::ref(LLG), std::ref(demag), calcDemag);
				writeLogStream(logstream, elapsed_ps, maximumTorque);
				nextLogOutput += prog.logStride;
			}
			outputTimer.add();
//////
			Map<MatrixXd_CM>(mag_vec.data(), nx, 3) = mag;
			odeTimer.start();
			// switching OFF ODEINT: Choice between CVODE and ODEINT is too complex, and ODEINT is incompatible with THRUST ///
			prog.useCVODE = true;
			if (prog.useCVODE) {
#ifdef USE_CUDA
				if (sd.useGPU && integrateOnGPU) {
					total_its += LLG.gpuODE(mag_vec, ode_start_t, ode_end_t, dt);
				} else {
#endif
					total_its += LLG.integrateSUNDIALS( mag_vec, ode_start_t, ode_end_t, dt );
#ifdef USE_CUDA
				}
#endif
			} else {
///// REMOVING GPU INTEGRATION WITH ODEINT (INCOMPATIBILITY THRUST / BOOST) ////
/*
#ifdef USE_CUDA 
				if (sd.useGPU && integrateOnGPU) {
					total_its += LLG.IntegrateOnGPU(mag_vec, ode_start_t, ode_end_t, dt);
				} else {
#endif
*/
					total_its += LLG.IntegrateODEINT(mag_vec, ode_start_t, ode_end_t, dt);
//#ifdef USE_CUDA 
////// 
//				}
//#endif
			}
			odeTimer.add();
			mag = Map<MatrixXd_CM>(mag_vec.data(), nx, 3);
			normalizeMag(mag);

			elapsed_t += ode_end_t;
			ode_start_t = elapsed_t;
			elapsed_ps = elapsed_t / realTimeScale;
			maximumTorque = LLG.getMaxTorque(mag);
			simTimer.end();
			double simRate = simTimer.getRate(ode_end_t / realTimeScale);
			outputTimer.start();
			std::cout << "elapsed time [ps]: " << std::fixed << std::setprecision(3) << elapsed_ps
//				<< ", iterations = " << total_its
					<< ", sim.rate [fs/s] = " << simRate << ", max.torque: " << std::scientific << maximumTorque << std::fixed
					<< "      \r" << std::flush;
			outputTimer.add();
			simTimer.start();
//		cudaProfilerStop();
			 totalTimer.add();
		}

//////////////////////////////////////////End of main loop /////////////////////////////////////////////////////////////

		// writing final output:
		calculateEnergies(std::ref(LLG), std::ref(demag), calcDemag);
		writeLogStream(logstream, elapsed_ps, maximumTorque);
		displayEnergies(elapsed_ps);
		if (prog.writerType == "GMV") {
			write.outputGMV(msh.name + ".gmr", mag);
		} else {
		write.outputVTK(msh.name, mag, demag.calcField(mag), "Magnetization", "Demag"); 
		}

		if (sd.Hys) {
			std::string hysOutfile = "hys_" + hys.getHysFileName(msh.name);
			write.outputVTK(hysOutfile, mag, sd.H_Hys, "external field");
			hys.moveFile(hysOutfile + ".vtu");
			hys.writeHysData(msh.NodeVolume, totalVolume, mag);
			writeLogStream(hys.hystLogFile, elapsed_ps, maximumTorque);
			std::cout << "HYSTERESIS: Finished calculation at "	<< sd.H_Hys * 1000. << " mT." << std::endl;
			sd.H_Hys = hys.nextField() / 1000.; // conversion mT to T
			if (hys.onBranch()) {
				LLG.setZeemanField(nx, sd.Hext, sd.theta_H, sd.phi_H, sd.H_Hys,	sd.theta_Hys, sd.phi_Hys);
				maximumTorque = LLG.getMaxTorque(mag);
				elapsed_t = elapsed_ps = 0.;
				nextGraphicsOutput = 0.;
				nextConsoleOutput = 0.;
				nextLogOutput = 0.;
			}
		}
	} while (sd.Hys && hys.onBranch());
	if (sd.Hys)
		hys.close();
	if (timeReached(elapsed_ps, (prog.duration + prog.initialTime))) {
		logstream << "# Maximum simulation time reached after " << std::fixed
				<< std::setprecision(3) << elapsed_ps << " picoseconds.\n";
	}
	if (maximumTorque <= prog.maxTorque) {
		logstream << "# Convergence reached: maximum torque " << maximumTorque
				<< " is below user-defined threshold.\n";
	}

	if (prog.showTimer) {
		showTimer();
		demag.outputTimer();
//		LLG.outputTimer();
	}
	writeLogFileFooter(logstream);
	logstream.close();
	std::cout << "Simulation finished." << std::endl;
//	std::cout << "Number of iterations: " << total_its << std::endl;

/////////// cleanup
	std::remove((msh.name + ".pid").c_str());
}

void TheSimulation::writeLogStream(std::ofstream &logstream, double elapsed_ps, double maximumTorque) {
	logstream <<  std::fixed << std::setprecision(4) << elapsed_ps << "\t" << totalEnergy << "\t" << DemagEnergy << "\t" << ExchangeEnergy << "\t"
	<< UniaxialEnergy << "\t" << ZeemanEnergy << "\t" << CubicEnergy << "\t" << SurfaceEnergy << "\t"
	<< DMIEnergy << "\t" << std::setprecision(3) << std::scientific << maximumTorque << "\t" << std::setprecision(10) << std::defaultfloat << getMeanM().transpose() <<"\t"
	<< std::fixed << std::setprecision(6) << meanH.transpose() * PhysicalConstants::mu0 << std::endl;
}

void TheSimulation::displayEnergies(double real_t) {
	std::cout << "                                                         \r";
	std::cout << "t[ps]: " << std::fixed << std::setprecision(3) << real_t << "\t"; // << std::setprecision(defaultPrecision);
	std::cout << "tot: " << totalEnergy << " de: " << DemagEnergy << " xc: "
			<< (std::abs(ExchangeEnergy) < 1.e-7 ? 0 : ExchangeEnergy)
			<< " un: " << UniaxialEnergy << " ze: "
			<< ZeemanEnergy << " cb: " << CubicEnergy << " sf: " << SurfaceEnergy << " dm: " << DMIEnergy 			
			<< "\tM: " << std::setprecision(5) << getMeanM().transpose() 
			<< std::endl; 
}

void TheSimulation::calculateEnergies( TheLLG& LLG, DemagField& demag, bool calcDemag ) {
	if (calcDemag)	 LLG.setHdem(demag.calcField(mag));
	DemagEnergy    = demag.getDemagEnergy(mag) / totalVolume;
	ExchangeEnergy = LLG.getExchangeEnergy(mag) / totalVolume;
	UniaxialEnergy = LLG.getUniaxialAnisotropyEnergy(mag) / totalVolume;
	ZeemanEnergy   = LLG.getZeemanEnergy(mag) / totalVolume;
	CubicEnergy	   = LLG.getCubicAnisotropyEnergy(mag) / totalVolume;
	SurfaceEnergy  = LLG.getSurfaceAnisotropyEnergy(mag) / totalVolume;
	DMIEnergy 	   = LLG.getDMIEnergy(mag) / totalVolume;
	totalEnergy = DemagEnergy + ExchangeEnergy + UniaxialEnergy + ZeemanEnergy + CubicEnergy + SurfaceEnergy + DMIEnergy;
	meanH 		=	LLG.getMeanH();	
}

void TheSimulation::generateInitialConfiguration() {
	InitialConfig initialState(msh.xyz, prog.initialConfiguration);
	mag = initialState.GetConfiguration();
	prog.resuming ? prog.initialTime = initialState.GetInitialTime() : prog.initialTime = 0.;
}

void TheSimulation::normalizeMag( MatrixXd& mag ) {
	if (sd.useGPU) {
#ifdef USE_CUDA 
		onGPU->NormalizeMag( mag, mag.rows() );
#endif
	} else {
		mag.rowwise().normalize();
	}
}

void TheSimulation::showTimer() {
	std::cout << "Timer output:\n";
	std::cout << "ode time [s]:\t" << odeTimer.durationInMus() / 1e6 << "\t(" << odeTimer.durationInMus() / totalTimer.durationInMus() * 100. <<" %)"<< std::endl;
	std::cout << "demag time [s]:\t" << pdeTimer.durationInMus() / 1e6 << "\t(" << pdeTimer.durationInMus() / totalTimer.durationInMus() * 100. <<" %)" << std::endl;
	std::cout << "output time [s]: " << outputTimer.durationInMus() / 1e6 << "\t(" << outputTimer.durationInMus() / totalTimer.durationInMus() * 100. <<" %)"<< std::endl;
	std::cout << "total time [s]:\t" << totalTimer.durationInMus() / 1e6  << std::endl;
}
