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
 * ProgramSpecs.cpp
 *
 *  Created on: May 3, 2017
 *      Author: riccardo
 */

#include "ProgramSpecs.h"
#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>
#include <iostream>
#include <fstream>
#include <tuple>
#include "auxiliaries.h"
#include <limits>
#include "PhysicalConstants.h"
#include "typedefs.h"

double inf = std::numeric_limits<double>::infinity();

enum Coords { x, y, z };

ProgramSpecs::ProgramSpecs() {
	configurationFile = "simulation.cfg";
	if (!fileExists(configurationFile)) noFileError();
	hasRequiredValues = true;
	initialize();
	firstFileNumber = 0;
}


void ProgramSpecs::initialize() {
	namespace opt = boost::program_options;
	opt::options_description desc("Options used in simulation.cfg:");
	desc.add_options()
			("use H2", opt::value<bool>()->default_value(true),"Use H2 matrix compression (logical)")
			("scale", opt::value<double>()->required(),	"Scaling factor relating a length unit to a distance in [m] (required)")
			("alpha", opt::value<double>()->required(),	"Gilbert damping constant (required).")
			("name", opt::value<std::string>()->required(), "Name of the project.")
			("initial state", opt::value<std::string>(), "initial configuration")
			("mesh type", opt::value<std::string>()->default_value("msh"), "mesh type")
			("resume", opt::value<bool>()->default_value(false), "continue a previous simulation")
			("torque limit", opt::value<double>()->default_value(1.e-4), "termination criterion: stop when torque is smaller than this")
			("first file number", opt::value<int>()->default_value(0), "number of first graphics output file")
			("solver type", opt::value<std::string>()->default_value("LU"), "select solver type: LU or CG")
			("CG tolerance", opt::value<double>()->default_value(1.e-7), "relative error tolerated in CG minimization")
			("remove precession", opt::value<bool>()->default_value(false), "Suppress precession, use only relaxation term of LLG")
			("duration", opt::value<double>()->default_value(inf), "Time after which the simulation ends")
			("time step", opt::value<double>()->default_value(0.5), "time step used in LLG integration")
			("log stride", opt::value<double>()->default_value(2), "output in log file")
			("config stride", opt::value<double>()->default_value(10), "output of spin configuration")
			("console stride", opt::value<double>()->default_value(10), "output in console")
			("writer type", opt::value<std::string>()->default_value("VTK"),"graphics output format: VTK or GMV.")
			("external field", opt::value<double>()->default_value(0.), "external field in mT.")
			("theta_H", opt::value<double>()->default_value(0.), "polar angle Hext direction (in degrees).")
			("phi_H", opt::value<double>()->default_value(0.), "azimuth angle Hext direction (in degrees).")
			("freeze demag", opt::value<bool>()->default_value(true), "separate demag field calculation from LLG integration (faster)")
			("neglect demag", opt::value<bool>()->default_value(false), "simulate without demag field")
			("remove demag", opt::value<bool>()->default_value(false), "simulate without demag field (alternative keyword)")
			("field pulse", opt::value<bool>()->default_value(false), "external field pulse applied (yes / no)")
			("pulse duration", opt::value<double>()->default_value(0.), "time width of Gaussian field pulse in ps")
			("pulse delay", opt::value<double>()->default_value(0.), "time in ps at which Gaussian field pulse reaches maximum value")
			("pulse amplitude", opt::value<double>()->default_value(0.), "peak value of Gaussian field pulse in mT")
			("pulse theta", opt::value<double>()->default_value(0.), "polar angle of field pulse direction (in degrees)")
			("pulse phi", opt::value<double>()->default_value(0.), "azimuthal angle of field pulse direction (in degrees)")
			("sweep field", opt::value<bool>()->default_value(false), "uses sweeping external field, changing linearly in time (yes / no)")
			("sweep start", opt::value<double>()->default_value(0.), "initial value of external field sweep in mT")
			("sweep end", opt::value<double>()->default_value(0.), "final value of external field sweep in mT")
			("sweep duration", opt::value<double>()->default_value(1.), "duration of the field sweeping in ps.")
			("sweep theta", opt::value<double>()->default_value(0.), "polar angle of sweeping field direction (in degrees)")
			("sweep phi", opt::value<double>()->default_value(0.), "azimuthal angle of sweeping field direction (in degrees)")
			("spin polarization", opt::value<double>()->default_value(0.), "degree of spin-polarization in stt dynamics")
			("current density", opt::value<double>()->default_value(0.), "current density in stt dynamics in [10^12 A/m^2]")
			("beta", opt::value<double>()->default_value(0.), "non-adiabaticity parameter in STT term")
			("current theta", opt::value<double>()->default_value(0.), "polar angle of electric current flow direction (in degrees)")
			("current phi", opt::value<double>()->default_value(0.), "azimuthal angle of electric current flow direction (in degrees)")
			("stt dynamics", opt::value<bool>()->default_value(false), "use STT term in dynamic calculation [yes /no]")
			("current pulse", opt::value<bool>()->default_value(false), "is the current applied as a pulse?")
			("current pulse duration", opt::value<double>()->default_value(1.), "time width of Gaussian current pulse in ps")
			("current pulse delay", opt::value<double>()->default_value(0.), "time in ps at which Gaussian current pulse reaches maximum value")
			("integrate on GPU",  opt::value<bool>()->default_value(true), "perform LLG integration on GPU instead of parallel CPU [yes /no]")
			("preconditioner", opt::value<std::string>()->default_value("ILU"), "preconditioner for GPU CG solver: ILU, MCILU, GS, MCGS")
			("current type", opt::value<std::string>()->default_value("none"), "options: PULSE, DC, or NONE")
			("local field", opt::value<std::string>()->default_value("none"), "application of local field: NONE, PULSE, SINE, CONST")
			("local field amplitude", opt::value<double>()->default_value(0), "maximum value of local field in mT")
			("local field frequency", opt::value<double>()->default_value(0), "frequency of local sinusoidal field in GHz")
			("local pulse duration", opt::value<double>()->default_value(1.), "time width of local field pulse in ps")
			("local pulse delay", opt::value<double>()->default_value(0), "time in ps at which local Gaussian field pulse is at maximum")
			("hysteresis", opt::value<bool>()->default_value(false), "calculate hysteresis through a sequence of states (yes / no)")
			("initial field", opt::value<double>()->default_value(0.), "initial field value of hysteresis [mT]")
			("final field", opt::value<double>()->default_value(0.), "final field value of hysteresis [mT]")
			("field step", opt::value<double>()->default_value(0.), "increment or decrement size in hysteresis [mT]")
			("hys theta", opt::value<double>()->default_value(0.), "polar angle of hysteresis field direction [degs]")
			("hys phi", opt::value<double>()->default_value(0.), "azimuthal angle of hysteresis field direction [degs]")
			("device number", opt::value<int>()->default_value(0), "identifier of GPU device in the case of multiple GPUs")
			("use CVODE", opt::value<bool>()->default_value(false), "select CVODE integrator (instead of ODEINT)")
			("timer output", opt::value<bool>()->default_value(false), "display amount of time spent in individual parts of the code")
			("gamma", opt::value<double>()->default_value(PhysicalConstants::gamma0), "gyromagnetic")
			;

	std::ifstream cfg_stream(configurationFile.c_str());
	try {
		bool allow_unregistered = true; // default boost value is false.
		opt::store(opt::parse_config_file(cfg_stream, desc, allow_unregistered), vm);
		notify(vm);
	} catch (const opt::error &ex) {
		 	std::cout << "ERROR - Input data in simulation.cfg: ";
		 	std::cerr << ex.what() << std::endl;
			exit(1);
	}
}


bool ProgramSpecs::noMissingValues(){return hasRequiredValues;}


void ProgramSpecs::evalInput() {
	if (defaultedCurrentType) {
		if (useSTT && !pulsedSTT) {
			std::cout << "simulation with DC current" << std::endl;
		}
		if (pulsedSTT) {
			std::cout << "simulation with pulsed current" << std::endl;
			useSTT = true;
		}
		return;
	}

	if (currentType == "dc") {
		if ( pulsedSTT ) {
			std::cout << "Conflicting information in simulation.cfg file:" << std::endl;
			std::cout << "current pulse = " << std::boolalpha << pulsedSTT << std::endl;
			std::cout << "current type = " << currentType << std::endl;
			exit(0);
		}
		std::cout << "simulation with DC current." << std::endl;
		useSTT = true;
		pulsedSTT = false;
	} else if (currentType == "pulse") {
		std::cout << "simulation with pulsed current" << std::endl;
		useSTT = true;
		pulsedSTT = true;
	} else if (currentType == "none") {
		if (pulsedSTT) {
			std::cout << "Conflicting information in simulation.cfg file:" << std::endl;
			std::cout << "current pulse = " << std::boolalpha << pulsedSTT << std::endl;
			std::cout << "current type = " << currentType << std::endl;
			exit(0);
		}
		if (useSTT) {
			std::cout << "Conflicting information in simulation.cfg file:" << std::endl;
			std::cout << "stt dynamics = " << std::boolalpha << useSTT << std::endl;
			std::cout << "current type = " << currentType << std::endl;
			exit(0);
				}
		useSTT = false;
		pulsedSTT = false;
	} else {
		std::cout << "Selected CURRENT TYPE (" << currentType << ") unknown." << std::endl;
		exit(0);
	}
}


void ProgramSpecs::readFile() {
	try { // will throw if required values are missing
		notify(vm);
	} catch (boost::program_options::error &e) {
		std::cout << "Error " << e.what() << std::endl;
		hasRequiredValues = false;
	}
	if (vm.count("scale")) scale = vm["scale"].as<double>();
	if (vm.count("use H2")) useH2 = vm["use H2"].as<bool>();
	if (vm.count("name")) name = vm["name"].as<std::string>();
	if (vm.count("alpha")) alpha = vm["alpha"].as<double>();
	if (vm.count("initial state")) initialConfiguration = vm["initial state"].as<std::string>();
	if (vm.count("resume")) resuming = vm["resume"].as<bool>();
	if (vm.count("torque limit")) maxTorque = vm["torque limit"].as<double>();
	if (vm.count("solver type")) solverType = toLower(vm["solver type"].as<std::string>());
	if (vm.count("CG tolerance")) cgTol = vm["CG tolerance"].as<double>();
	if (vm.count("remove precession")) noPrecession = vm["remove precession"].as<bool>();
	if (vm.count("duration")) duration = vm["duration"].as<double>();
	if (vm.count("time step")) timeStep = vm["time step"].as<double>();
	if (vm.count("log stride")) logStride = vm["log stride"].as<double>();
	if (vm.count("config stride")) cfgStride = vm["config stride"].as<double>();
	if (vm.count("console stride")) consoleStride = vm["console stride"].as<double>();
	if (vm.count("writer type")) writerType = toLower(vm["writer type"].as<std::string>());
	if (vm.count("external field")) Hext = vm["external field"].as<double>();
	if (vm.count("theta_H")) theta_H = vm["theta_H"].as<double>();
	if (vm.count("phi_H")) phi_H = vm["phi_H"].as<double>();
	if (vm.count("freeze demag")) freezeDemag = vm["freeze demag"].as<bool>();
	if (vm.count("mesh type")) meshType = toLower(vm["mesh type"].as<std::string>());
	if (vm.count("field pulse")) fieldPulse = vm["field pulse"].as<bool>();
	if (vm.count("pulse duration")) pulseWidth = vm["pulse duration"].as<double>();
	if (vm.count("pulse delay")) pulseDelay = vm["pulse delay"].as<double>();
	if (vm.count("pulse amplitude")) pulsePeak = vm["pulse amplitude"].as<double>();
	if (vm.count("pulse theta")) pulseTheta = vm["pulse theta"].as<double>();
	if (vm.count("pulse phi")) pulsePhi = vm["pulse phi"].as<double>();
	if (vm.count("remove demag")) noDemag = (vm["neglect demag"].as<bool>() || vm["remove demag"].as<bool>());
	if (vm.count("sweep field")) sweepUsed = vm["sweep field"].as<bool>();
	if (vm.count("sweep start")) sweepStart = vm["sweep start"].as<double>();
	if (vm.count("sweep end")) sweepEnd = vm["sweep end"].as<double>();
	if (vm.count("sweep duration")) sweepDuration = vm["sweep duration"].as<double>();
	if (vm.count("sweep theta")) sweepTheta = vm["sweep theta"].as<double>();
	if (vm.count("sweep phi")) sweepPhi = vm["sweep phi"].as<double>();
	if (vm.count("spin polarization")) P = vm["spin polarization"].as<double>();
	if (vm.count("current density")) jc_val = vm["current density"].as<double>();
	if (vm.count("beta")) beta = vm["beta"].as<double>();
	if (vm.count("current theta")) theta_j = vm["current theta"].as<double>();
	if (vm.count("current phi")) phi_j = vm["current phi"].as<double>();
	if (vm.count("stt dynamics")) useSTT = vm["stt dynamics"].as<bool>();
	if (vm.count("integrate on GPU")) integrateOnGPU = vm["integrate on GPU"].as<bool>();
	if (vm.count("preconditioner")) preconditionerType = toLower(vm["preconditioner"].as<std::string>());
	if (vm.count("current pulse")) pulsedSTT = vm["current pulse"].as<bool>();
	if (vm.count("current pulse delay")) sttPulseDelay = vm["current pulse delay"].as<double>();
	if (vm.count("current pulse duration")) sttPulseWidth = vm["current pulse duration"].as<double>();
	if (vm.count("current type")) currentType = toLower(vm["current type"].as<std::string>());
	if (vm.count("local field")) localFieldType = toLower(vm["local field"].as<std::string>());
	if (vm.count("local field amplitude")) localAmplitude = vm["local field amplitude"].as<double>();
	if (vm.count("local field frequency")) localFrequency = vm["local field frequency"].as<double>();
	if (vm.count("local pulse duration")) localPulseWidth = vm["local pulse duration"].as<double>();
	if (vm.count("local pulse delay")) localPulseDelay = vm["local pulse delay"].as<double>();
	if (vm.count("hysteresis")) useHyst = vm["hysteresis"].as<bool>();
	if (vm.count("initial field")) firstField = vm["initial field"].as<double>();
	if (vm.count("final field")) lastField = vm["final field"].as<double>();
	if (vm.count("hys theta")) hystTheta = vm["hys theta"].as<double>();
	if (vm.count("hys phi")) hystPhi = vm["hys phi"].as<double>();
	if (vm.count("field step")) deltaH = vm["field step"].as<double>();
	if (vm.count("device number")) deviceNumber = vm["device number"].as<int>();
	if (vm.count("use CVODE"))useCVODE = vm["use CVODE"].as<bool>();
	if (vm.count("timer output"))showTimer = vm["timer output"].as<bool>();
	if (vm.count("gamma"))gamma=vm["gamma"].as<double>();

	defaultedCurrentType = vm["current type"].defaulted();
	evalInput();
}


std::string ProgramSpecs::getName() {return name;}


double ProgramSpecs::getExternalField() {return Hext;}


std::tuple<double, double> ProgramSpecs::getHextAngles(){
	return std::make_tuple(theta_H, phi_H);
}


void ProgramSpecs::noFileError() {
	std::cout << "Configuration file 'simulation.cfg' not found.\n" << std::endl;
	exit(0);
}


void Hlocal::setProfile(const Eigen::MatrixXd& Profile_) {
	Profile = Profile_;
}


Eigen::MatrixXd Hlocal::gaussPulseValue(const double t) {
	const double sigma = (t - pulseDelay) / pulseWidth;
	const double xpnt = (sigma * sigma) / 2.;
	return Profile * amplitude * std::exp(-xpnt);
} 


Eigen::MatrixXd Hlocal::sineFieldValue(const double t) {
	return Profile * amplitude * std::cos(omega * t);
}


Eigen::MatrixXd Hlocal::constantValue(const double) {
	return Profile * amplitude;
}


void Hlocal::setValueFunction() {
	if (type == "pulse") {
		localField = [this](const double t)-> Eigen::MatrixXd {return gaussPulseValue(t);};
	} else if (type == "sine") {
		localField = [this](const double t)-> Eigen::MatrixXd {return sineFieldValue(t);};
	} else {
		localField = [this](const double t)-> Eigen::MatrixXd {return constantValue(t);};
	}
}


void Hlocal::printData() {
	std::cout << "type: " << type << "\n";
	std::cout << "amplitude: " << amplitude * PhysicalConstants::mu0 * 1000. << " mT\n";
	std::cout << "frequency: " << frequency << " ps^-1" << std::endl;
	std::cout << "delay: " << pulseDelay << " ps" << std::endl;
	std::cout << "duration: " << pulseWidth << " ps" << std::endl;
	std::cout << "use local data: " << std::boolalpha << isUsed << std::endl;
	std::cout << "Profile, row 42: " << Profile.row(42) << std::endl;
}


Eigen::MatrixXd Hlocal::getSpatialProfile() {
	return Profile;
}


Hlocal ProgramSpecs::getHloc() {
	hl.amplitude = localAmplitude;
	hl.frequency = localFrequency;
	
	hl.pulseDelay = localPulseDelay;
	hl.pulseWidth = localPulseWidth;
	hl.type = localFieldType;
	hl.isUsed = (hl.type == "const" || hl.type == "pulse" || hl.type == "sine");
	if (!hl.isUsed && hl.type != "none") {
		std::cout << "Local field option '" << hl.type << "' not understood." << std::endl;
		std::cout << "Not using any local field." << std::endl;
	}
	hl.getSpatialProfile();
	return hl;
}


Hdynamic ProgramSpecs::getHdyn() {
//	enum Coords {	x, y, z 	};
	Hdynamic hp;
	hp.pulseDelay = pulseDelay;
	hp.pulsePeak = pulsePeak;
	hp.width = pulseWidth;
	hp.pulseIsUsed = fieldPulse;
	hp.pulseTheta = pulseTheta * PhysicalConstants::pi / 180.;
	hp.pulsePhi = pulsePhi * PhysicalConstants::pi / 180.;
	hp.sweepIsUsed = sweepUsed;
	hp.sweepStart = sweepStart;
	hp.sweepEnd = sweepEnd;
	hp.sweepDuration = sweepDuration;
	hp.sweepTheta = sweepTheta * PhysicalConstants::pi / 180.;
	hp.sweepPhi = sweepPhi * PhysicalConstants::pi / 180.;

/*
	hp.pulseDir(x)= std::sin(hp.theta) * std::cos(hp.phi);
	hp.pulseDir(y)= std::sin(hp.theta) * std::sin(hp.phi);
	hp.pulseDir(z)= std::cos(hp.theta);
	*/
	return hp;
}


double Hdynamic::gaussPulseValue(double t) {
		const double sigma = ( t - pulseDelay ) / width;
		const double xpnt = (sigma * sigma) / 2.;
		return pulsePeak * std::exp(-xpnt);
}


double Hdynamic::sweepFieldValue(double t) {
		double s = sweepStart;
		double e = sweepEnd;
		double d = sweepDuration;
		double slope = ( e - s ) / d ;
		double fieldValue = sweepStart + slope * t;
		if ((fieldValue > sweepEnd  && slope > 0.) || (fieldValue < sweepEnd  && slope < 0.)) fieldValue = sweepEnd;
		return fieldValue;
	}


STT::STT(SpMat gradX_, SpMat gradY_, SpMat gradZ_) :
			gradX(gradX_), gradY(gradY_), gradZ(gradZ_){}


void STT::setValues(int nx, double scale, const ProgramSpecs& prog) {

	double P = prog.P;
	double jc_val = prog.jc_val;
	beta = prog.beta; // also known as xi
	double theta = prog.theta_j;
	double phi = prog.phi_j;

	pulseIsUsed = prog.pulsedSTT;
	pulseWidth = prog.sttPulseWidth;
	pulseDelay = prog.sttPulseDelay;

	jc_vec = Eigen::MatrixXd::Ones(nx,3);
	double deg2rad = PhysicalConstants::pi / 180.;
	jc_vec.col(x) = Eigen::VectorXd::Ones(nx) * std::sin(theta * deg2rad) * std::cos(phi*deg2rad);
	jc_vec.col(y) = Eigen::VectorXd::Ones(nx) * std::sin(theta * deg2rad) * std::sin(phi*deg2rad);
	jc_vec.col(z) = Eigen::VectorXd::Ones(nx) * std::cos(theta * deg2rad);
	jc_vec *= (1.e12 * jc_val);
	btimesMs = P * PhysicalConstants::mub  * PhysicalConstants::mu0 / (PhysicalConstants::gamma0 * PhysicalConstants::e_charge * ( 1 + beta * beta )) / scale;
	Ustt.resize(nx, 3);
}


void STT::setEta(const Eigen::VectorXd& invJs) {
	Eigen::VectorXd eta = btimesMs * invJs;
	eta_jx = jc_vec.col(x).cwiseProduct(eta);
	eta_jy = jc_vec.col(y).cwiseProduct(eta);
	eta_jz = jc_vec.col(z).cwiseProduct(eta);
};


double STT::gaussPulseValue(double t) {
	double s = pulseWidth;
	double a = 1.; // pulsePeak;
	double b = pulseDelay;
	double xpnt = ((t-b)/s * (t-b)/s) / 2.;
	double p = a * std::exp(-xpnt);
	return p;
	return 0;
}
