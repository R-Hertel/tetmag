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
 * DemagField.h
 *
 *  Created on: May 4, 2017
 *      Author: riccardo
 */

#ifndef DEMAGFIELD_H_
#define DEMAGFIELD_H_
#include <Eigen/SparseCore>
#include <Eigen/Dense>
#include <vector>
#include <string>
#include <functional> // for std::function
#include "SolverFactory.h"
#include "MeshData.h"
#ifdef USE_CUDA
#include "GpuPotential.h"
#endif
#include "typedefs.h"
#include "Timer.h"
#include <memory> // for std:::shared_ptr

typedef double* pvector;

class DemagField {
private:
	size_t nx;
	size_t bnx;
	std::vector<int> boundaryNodes;
	Eigen::MatrixXd dirichletBEM;
	SpMat dirichletMatrix;
	SpMat neumannMatrix;
	int fixedNeumannNode; // value must be set at one node since the solution of the pure Neumann problem is not unique
	SpMat tGradX, tGradY, tGradZ;
	SpMat gradX, gradY, gradZ;
	Eigen::VectorXd Js;
	Eigen::VectorXd NodeVolume;
	bool useH2;
	Eigen::MatrixXd Hdem;
	void prepareNeumannMatrix();
	Eigen::VectorXd rhsPoisson(MRef &);
	Eigen::VectorXd solvePoisson(const Eigen::Ref<const Eigen::VectorXd>&);
	Eigen::VectorXd solveLaplace(const Eigen::Ref<const Eigen::VectorXd>&);
	Eigen::VectorXd boundaryIntegral(const Eigen::Ref<const Eigen::VectorXd>&);
	Eigen::VectorXd h2MVP(Eigen::VectorXd&);
//	Eigen::VectorXd hMVP(Eigen::VectorXd&);
	Eigen::VectorXd denseMVP(Eigen::VectorXd&);
	std::function <Eigen::VectorXd ( Eigen::VectorXd& ) > selectedMVPtype;
	Eigen::VectorXd mvp(Eigen::VectorXd&);
	std::shared_ptr<SolverFactory> dirichletSolver;
	std::shared_ptr<SolverFactory> neumannSolver;
	Eigen::VectorXd calcPotential(MRef &);
	Eigen::MatrixXd calcGradientField(const Eigen::Ref<const Eigen::VectorXd>&);
	bool useGPU;
	double cgTol;

	Timer u1Timer, h2Timer, u2Timer, demagTimer;
#ifdef USE_CUDA
	GpuPotential gpuField;
#endif
public:
	void outputTimer();
	Eigen::VectorXd getDivM(MRef &);
	DemagField(const MeshData&, const Eigen::VectorXd& , bool);
	double getDemagEnergy(const Eigen::Ref<const Eigen::MatrixXd>&);
	Eigen::MatrixXd calcField(MRef &);
	void initializeSolvers(std::string, std::string, double);
	DemagField(){};
};
#endif /* DEMAGFIELD_H_ */
