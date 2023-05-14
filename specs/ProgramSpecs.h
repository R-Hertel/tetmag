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
 * ProgramSpecs.h
 *
 *  Created on: May 3, 2017
 *      Author: riccardo
 */

#ifndef PROGRAMSPECS_H_
#define PROGRAMSPECS_H_
#include <string>
#include <boost/program_options.hpp>
#include <tuple>
#include <Eigen/Dense>
#include "typedefs.h"
#include <functional>

struct Hdynamic {
	double width;
	double pulseDelay;
	double pulsePeak;
	bool pulseIsUsed;
	double gaussPulseValue(double);
	double pulseTheta;
	double pulsePhi;
	bool sweepIsUsed;
	double sweepTheta;
	double sweepPhi;
	double sweepStart;
	double sweepEnd;
	double sweepDuration;
	Eigen::MatrixXd pulseH;
	double sweepFieldValue(double);
	Eigen::MatrixXd sweepH;
};

class Hlocal {
private:
	Eigen::MatrixXd gaussPulseValue(const double);
	Eigen::MatrixXd sineFieldValue(const double);
	Eigen::MatrixXd constantValue(const double);
	Eigen::MatrixXd Profile;
public:
	double omega;
	bool isUsed;
	std::string type;
	double frequency;
	double amplitude;
	double pulseWidth;
	double pulseDelay;
	std::function <Eigen::MatrixXd(const double)> localField;
	void setValueFunction();
	void setProfile(const Eigen::MatrixXd&);
	void printData();
	Eigen::MatrixXd getSpatialProfile();
};

class ProgramSpecs;
//class SimulationData;

struct STT {
	double pulseDelay;
	double pulseWidth;
	bool pulseIsUsed;
	SpMat gradX, gradY, gradZ;
	Eigen::MatrixXd jc_vec;
	double beta; // also known as xi
	double btimesMs;
	Eigen::VectorXd eta;
	Eigen::VectorXd eta_jx, eta_jy, eta_jz;
	STT(){};
	STT(SpMat, SpMat, SpMat);
	Eigen::MatrixXd Ustt;
	void setValues(int nx, double scale, const ProgramSpecs&);
	void setEta(const Eigen::VectorXd&);
	double gaussPulseValue(double);
};

class ProgramSpecs {
private:
	std::string configurationFile;
	boost::program_options::variables_map vm;
	std::string name;
	void initialize();
	void noFileError();
	bool hasRequiredValues;
	double Hext;
	double theta_H;
	double phi_H;
	double pulseWidth;
	double pulseDelay;
	double pulsePeak;
	double pulseTheta;
	double pulsePhi;
	double sweepStart;
	double sweepEnd;
	double sweepTheta;
	double sweepPhi;
	double sweepDuration;
	void evalInput();
	double localPulseDelay;
	double localPulseWidth;
	double localFrequency;
	double localAmplitude;
	std::string localFieldType;
public:
	bool fieldPulse;
	bool sweepUsed;
	std::string currentType;
	bool defaultedCurrentType;
	bool pulsedSTT;
	double sttPulseDelay;
	double sttPulseWidth;
	double alpha;
	double scale;
	bool noMissingValues();
	std::string getName();
	double getExternalField();
	std::tuple<double, double> getHextAngles();
	ProgramSpecs();
	void readFile();
	std::string initialConfiguration;
	std::string solverType;
	std::string preconditionerType;
	int deviceNumber;
	std::string writerType;
	std::string localField;
	bool useH2;
	bool resuming;
	double initialTime;
	int firstFileNumber;
	bool noPrecession;
	double duration;
	double timeStep;
	double logStride;
	double cfgStride;
	double consoleStride;
	double maxTorque;
	bool freezeDemag;
	bool noDemag;
	std::string meshType;
	Hdynamic getHdyn();
	Hlocal getHloc();
	double P;
	double jc_val;
	double beta;
	double theta_j;
	double phi_j;
	bool useSTT;
	bool integrateOnGPU;
	bool useHyst;
	double firstField, lastField;
	double hystPhi, hystTheta;
	double hystStart, hystEnd, deltaH;
	Hlocal hl;
	double cgTol;
	bool useCVODE;
	bool showTimer;
};


#endif /* PROGRAMSPECS_H_ */
