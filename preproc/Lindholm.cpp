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
 * Lindholm.cpp
 *
 *  Created on: Feb 25, 2019
 *      Author: riccardo
 */

#include <Eigen/Dense>
#include "Lindholm.h"
#include "auxiliaries.h"
#include <cmath>
//#include "PhysicalConstants.h"
#include <iostream>
using namespace Eigen;


double Lindholm::solidAngle(const Vector3d& x, const Vector3d& T1, const Vector3d& T2, const Vector3d& T3) {
	Vector3d r1, r2, r3;
	r1 = T1 - x;
	r2 = T2 - x;
	r3 = T3 - x;
	double r1_n = r1.norm();
	double r2_n = r2.norm();
	double r3_n = r3.norm();
	double r123 = r1_n * r2_n * r3_n;
	double r23 = r2.dot(r3);
	double r12 = r1.dot(r2);
	double r31 = r3.dot(r1);
	double numer = r123 + r1_n * r23 + r2_n * r31 + r3_n * r12;
	double denom = std::sqrt(
			2.0 * (r2_n * r3_n + r23) * (r3_n * r1_n + r31)
					* (r1_n * r2_n + r12));
	double arg = numer/denom;
	if (areEqual( arg, 1.)) { arg = 1.; } // roundoff errors can lead to arguments > 1 resulting in NaN
	double omega = 2.0 * std::acos(double(arg));
	return omega;
}


Vector3d Lindholm::weights(const Vector3d& x0, const Vector3d& p1, const Vector3d& p2,
		const Vector3d& p3, double A, const Vector3d& nv, bool isOrdered) {
	// Implementation of D.A. Lindholm, IEEE Trans. Mag 20, 2025 (1984)
	MatrixXd p(3, 3);
	p.col(0) = p1;
	p.col(1) = p2;
	p.col(2) = p3;
	bool containsX0 = ((p1 - x0).isMuchSmallerThan(1.)
			|| (p2 - x0).isMuchSmallerThan(1.)
			|| (p3 - x0).isMuchSmallerThan(1.));
	// do not use .isApprox() here. That comparison is not reliable if x0 is the zero vector.
	if (containsX0) {
		return Vector3d::Zero();
	}

	Vector3d chi = isOrdered? nv : triangleUnitNormVector(p.col(0), p.col(1), p.col(2));
	MatrixXd rho_vec(3, 3);
	for (int i = 0; i < 3; ++i) {
		rho_vec.col(i) = p.col(i) - x0;
	}
	MatrixXd xi(3, 3);
	for (int i = 0; i < 3; ++i) {
		xi.col(i) = rho_vec.col((i + 1) % 3) - rho_vec.col(i);
	}
	VectorXd eta0(3);
	VectorXd s(3);
	MatrixXd eta(3, 3);
	VectorXd rho(3);

	for (int i = 0; i < 3; ++i) {
		s(i) = xi.col(i).norm();
		xi.col(i) = static_cast<Vector3d>(xi.col(i)).normalized();
		eta.col(i) = chi.cross(static_cast<Vector3d>(xi.col(i))); // this is already normalized
		eta0(i) = eta.col(i).dot(rho_vec.col(i));
		rho(i) = rho_vec.col(i).norm();
	}

	double chi0 = chi.dot(rho_vec.col(0));
	MatrixXd gamma(3, 3);
	for (int i = 0; i < 3; ++i) {
		for (int j = 0; j < 3; ++j) {
			gamma(i, j) = xi.col((i + 1) % 3).dot(xi.col(j));
		}
	}

	Vector3d P;
	for (int i = 0; i < 3; ++i) {
		P[i] = std::log(
				(rho[i] + rho[(i + 1) % 3] + s[i])
						/ (rho[i] + rho[(i + 1) % 3] - s[i]));
	}

	double Omega = solidAngle(x0, p.col(0), p.col(1), p.col(2));
	if (chi0 < 0.)
		Omega = -Omega;
	Vector3d L = Vector3d::Zero();
	for (int i = 0; i < 3; ++i) {
		L[i] = s[(i + 1) % 3] / (8. * pi * A)
				* (eta0[(i + 1) % 3] * Omega - chi0 * (gamma * P)[i]);
	}
	if (!isOrdered && nv.dot(chi) < 0) {
		L = -L;
	}
	return L;
}
