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
#include "SolverFactory.h"
#include "MeshData.h"
#include "BEMOperator.h"
#include "DemagBackend.h"
#include "typedefs.h"
#include "Timer.h"
#include <memory> // for std:::shared_ptr

class DemagField {
private:
	size_t nx;
	SpMat dirichletMatrix;
	SpMat neumannMatrix;
	int fixedNeumannNode; // value must be set at one node since the solution of the pure Neumann problem is not unique
	SpMat tGradX, tGradY, tGradZ;
	SpMat gradX, gradY, gradZ;
	Eigen::VectorXd Js;
	Eigen::VectorXd NodeVolume;
	void prepareNeumannMatrix();
	std::shared_ptr<SolverFactory> makeAndSetupSolver(const std::string&, const std::string&, const SpMat&, const std::string&);
	std::shared_ptr<BEMOperator> bem;
	std::shared_ptr<DemagBackend> backend;
	double cgTol;

	Timer u1Timer, h2Timer, u2Timer, demagTimer;
public:
	void outputTimer();
	DemagField(const MeshData&, const Eigen::VectorXd& , bool);
	double getDemagEnergy(const Eigen::Ref<const Eigen::MatrixXd>&);
	Eigen::MatrixXd calcField(MRef &);
	void initializeSolvers(std::string, std::string, double);
	DemagField(){};
};
#endif /* DEMAGFIELD_H_ */
