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
 * Materials.cpp
 *
 *  Created on: Mar 21, 2017
 *      Author: riccardo
 */
#include "Materials.h"
#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <cmath>
#include "PhysicalConstants.h"
using namespace Eigen;

Material::Material(){}


Material::Material(std::string desc_) {
	if (desc_ == "Permalloy") {
		A = 1.3e-11; Js = 1.0;	Ku = 500.; phi_u = 90.; theta_u = 0.; Kc1 = 0.; Kc2 = 0.;
		phi_Euler = 0.;	theta_Euler = 0.; psi_Euler = 0.; Ks = 0.;
		desc = desc_;
		D = 0.;
	} else if (desc_ == "Cobalt") {
		A = 2.5e-11; Js = 2.1; Ku = 3000.; phi_u = 90.; theta_u = 0.; Kc1 = 0.; Kc2 = 0.;
		phi_Euler = 0.;	theta_Euler = 0.; psi_Euler = 0.; Ks = 3.8e-3;
		desc = desc_;
		D = 0.;
	} else {
		desc = "Iron";
		A = 2.5e-11; Js = 1.6;	Ku = 1500.; phi_u = 90.; theta_u = 0.; Kc1 = 1000.; Kc2 = 0.;
		phi_Euler = 0.; theta_Euler = 0.; psi_Euler = 0.; Ks = 0.;
		D = 0.;
	}
}


Material::Material(std::string desc_, double A_, double Js_, double Ku_,
		double theta_u_, double phi_u_, double Kc1_, double Kc2_, double Ks_, double phi_c_,
		double theta_c_, double psi_c_, double D_ = 0.) :
				desc(desc_), A(A_), Js(Js_), Ku(Ku_), theta_u(theta_u_), phi_u(phi_u_),
		Kc1(Kc1_), Kc2(Kc2_), Ks(Ks_), phi_Euler(phi_c_), theta_Euler(theta_c_),
		psi_Euler(psi_c_), D(D_){}


Matrix3d Material::calculateCubicAxesVectors() {
	double pi = PhysicalConstants::pi;
	double gradToRad = pi / 180.;
	Matrix3d axes;
	axes =  AngleAxisd(phi_Euler * gradToRad, Vector3d::UnitZ());
	axes =	axes * AngleAxisd(theta_Euler * gradToRad, Vector3d::UnitX());
	axes =  axes * AngleAxisd( psi_Euler  * gradToRad, Vector3d::UnitZ());
//  The code corresponds to a chaining like AngleAxisd(alpha..) * AngleAxisd(beta..) * AngleAxisd(gamma..).
//  The order of the multiplication sequence is important since the rotations are not commutative.
//
//  The axes are stored in the COLS of the Matrix: col(0) is the rotated x-axis, col(1) is the y-axis, col(2) the z-axis.
//  Using the "x-convention" of the Euler angles: First rotate by phi around z,
//  then rotate by theta around x', and finally a rotation by psi around z'.
//
//  Example: the angles phi = 45°, theta = 54.74° rotate the z axis to the <111> direction
	cubicAxes = axes.unaryExpr([](double a) {double small = 1.e-14; return a * (std::abs(a) > small);});
	// This ^^^ sets ridiculously small numerical values to zero. Needs -std=c++11 compilation because of lambda.
	return cubicAxes;
}


VectorXi Material::assignMaterials(const MatrixXd& xyz) {
	enum mats {one, two, three};
	enum coords {x, y, z};
	int nx = xyz.rows();
	VectorXi NodeMaterial = VectorXi::Zero(nx);
/*
// position-dependent material assignment:
	for (int i = 0; i < nx; ++i) {
		if (xyz(i, z) < 0.1)
			NodeMaterial(i) = one;
		else
			NodeMaterial(i) = two;
	}
*/
	return NodeMaterial;
}


