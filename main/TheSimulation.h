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
 * TheSimulation.h
 *
 *  Created on: Nov 30, 2016
 *      Author: riccardo
 */

#ifndef THESIMULATION_H_
#define THESIMULATION_H_

#include "SimulationData.h"
#include "MeshData.h"
#include "ProgramSpecs.h"
#include "TheLLG.h"
#include "DemagField.h"
#include <fstream>
#include <string>
#include <Eigen/Dense>
#include <Eigen/SparseCore>
#include "typedefs.h"
#include <memory>
#include "Timer.h"

#ifdef USE_CUDA
	#include "EffFieldGPU.h"
//	#include <paralution/paralution.hpp>
	#include "LLGWrapper.h"
#endif

class TheSimulation {
private:
	size_t nx;
	SimulationData sd;
	MeshData msh;
	ProgramSpecs prog;
	void displayEnergies(double);
	void assignMaterialsToFE();
	Eigen::Vector3d getMeanM();
	Eigen::MatrixXd mag;
	double DemagEnergy;
	double ExchangeEnergy;
	double UniaxialEnergy;
	double ZeemanEnergy;
	double CubicEnergy;
	double SurfaceEnergy;
	double DMIEnergy;
	double totalEnergy;
	void calculateEnergies(TheLLG&, DemagField&, bool);
	double totalVolume;
/*	
	double directDMI;
	double directExch;
*/	
//	double energyScale;
	void writeLogStream(std::ofstream&, double, double);
	TheLLG AlternateLLG(bool, int, DemagField&);
	Eigen::Vector3d meanH;
	void normalizeMag(Eigen::MatrixXd&);
	Timer simTimer, odeTimer, pdeTimer, outputTimer, totalTimer;
	void showTimer();
#ifdef USE_CUDA
	std::shared_ptr<EffFieldGPU> onGPU;
	LLGWrapper gpuWrap;
#endif

public:
	TheSimulation(SimulationData&, MeshData&, ProgramSpecs&);
	void start();
	void generateInitialConfiguration();
};

#endif /* THESIMULATION_H_ */
