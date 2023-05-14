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
 * Materials.h
 *
 *  Created on: Mar 21, 2017
 *      Author: riccardo
 */

#ifndef MATERIALS_H_
#define MATERIALS_H_
#include <string>
#include <Eigen/Dense>
#include <Eigen/Geometry>

class Material {
public:
	std::string desc; // description of the material
	double A;  // exchange constant
	double Js; // saturation magnetic polarization Js = mu0 * Ms
	double Ku; // uniaxial anisotropy constant in J/m^3
	double theta_u;
	double phi_u; // direction of uniaxial easy axis, measured in degs
	double Kc1;   // first cubic anisotropy constant in J/m^3
	double Kc2;   // second cubic anisotropy constant
	double Ks; // surface anisotropy constant
	double phi_Euler; // Euler angles
	double theta_Euler;
	double psi_Euler;
	Eigen::Matrix3d cubicAxes; // orthonormal set of easy axes of cubic anisotropy
	Material(std::string, double, double, double, double, double, double, double, double, double, double, double, double);
	Material(std::string);
	Material();
	Eigen::Matrix3d calculateCubicAxesVectors();
	double D; // DMI constant in J/m2 (volume term)
	Eigen::VectorXi assignMaterials(const Eigen::MatrixXd&);
};

#endif /* MATERIALS_H_ */
