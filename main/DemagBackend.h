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
 * DemagBackend.h
 *
 *  Created on: Aug 7, 2026
 *      Author: riccardo
 */

#ifndef MAIN_DEMAGBACKEND_H_
#define MAIN_DEMAGBACKEND_H_
#include <Eigen/Dense>
#include <functional>
#include "typedefs.h"

// Interface to a demagnetizing-field implementation. 
// The solver stages carry no state in their signatures. 
// Each backend owns the magnetization and the resulting field.

struct DemagBackend
{
    std::function<void(MRef &)> setMagnetization;
    std::function<void()> computeRhs;
    std::function<void()> solvePoisson;
    std::function<const Eigen::VectorXd &()> poissonPotential;
    std::function<void(const Eigen::Ref<const Eigen::VectorXd> &)> setBoundaryValues;
    std::function<void()> solveLaplace;
    std::function<void()> computeField;
    std::function<const Eigen::MatrixXd &()> field;
};

#endif /* MAIN_DEMAGBACKEND_H_ */
