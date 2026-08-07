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
 * BEMOperator.h
 *
 *  Created on: Aug 7, 2026
 *      Author: riccardo
 */

#ifndef MAIN_BEMOPERATOR_H_
#define MAIN_BEMOPERATOR_H_
#include <Eigen/Dense>
#include <functional>
#include <vector>
#include "MeshData.h"

typedef double* pvector;

class BEMOperator {
private:
	size_t nx;
	size_t bnx;
	std::vector<int> boundaryNodes;
	Eigen::MatrixXd dirichletBEM;
	Eigen::VectorXd h2MVP(Eigen::VectorXd&);
	Eigen::VectorXd denseMVP(Eigen::VectorXd&);
	std::function<Eigen::VectorXd(Eigen::VectorXd&)> mvp;
public:
	BEMOperator(const MeshData&, size_t, bool);
	Eigen::VectorXd boundaryIntegral(const Eigen::Ref<const Eigen::VectorXd>&);
};

#endif /* MAIN_BEMOPERATOR_H_ */
