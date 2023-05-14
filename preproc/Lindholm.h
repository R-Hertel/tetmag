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
 * Lindholm.h
 *
 *  Created on: Feb 25, 2019
 *      Author: riccardo
 */

#ifndef LINDHOLM_H_
#define LINDHOLM_H_
#include <Eigen/Dense>
#include "PhysicalConstants.h"

class Lindholm {
private:
	double pi = PhysicalConstants::pi;
public:
	double solidAngle(const Eigen::Vector3d& x, const Eigen::Vector3d& T1, const Eigen::Vector3d& T2, const Eigen::Vector3d& T3);
	Eigen::Vector3d weights(const Eigen::Vector3d&, const Eigen::Vector3d&, const Eigen::Vector3d&, const Eigen::Vector3d&, double, const Eigen::Vector3d&, bool isOrdered = false);
};

#endif /* LINDHOLM_H_ */
