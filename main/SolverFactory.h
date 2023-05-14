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
 * SolverFactory.h
 *
 *  Created on: Mar 7, 2017
 *      Author: riccardo
 */

#ifndef SOLVERFACTORY_H_
#define SOLVERFACTORY_H_
#include <Eigen/SparseCore>
#include <string>
#include "typedefs.h"
#include <memory>

typedef Eigen::SparseMatrix<double, Eigen::RowMajor> SpMat_RM;
typedef Eigen::SparseMatrix<double, Eigen::ColMajor> SpMat_CM;

class SolverFactory {
protected:
	double cgTol;
	SpMat A;
	SpMat_RM Ar;
	SpMat_CM Ac;
	Eigen::VectorXd b;
	Eigen::VectorXd x;
	std::string type;
	Eigen::VectorXd previous;
	std::string preconditionerType;
public:
	static std::shared_ptr<SolverFactory> makeSolver(std::string);
	void setLoadVector(const Eigen::Ref<const Eigen::VectorXd>&);
	virtual void solve() = 0;
	virtual bool compute() = 0;
	virtual void report() = 0;
	virtual bool wasSuccessful() = 0;
	virtual ~SolverFactory() {}
	Eigen::VectorXd result();
	void setMatrix(const SpMat&  );
	void definePreconditioner(const std::string );
	void setCGTolerance(const double);
};

#endif /* SOLVERFACTORY_H_ */
