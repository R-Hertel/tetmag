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
 * InitMag.cpp
 *
 *  Created on: Oct 12, 2016
 *      Author: riccardo
 */
#include "InitMag.h"
#include "VTKReader.h"
#include "PhysicalConstants.h"
#include <cmath>
#include <random>
#include <Eigen/Dense>
#include <boost/algorithm/string.hpp>
#include <iostream>
#include <random>
#include "auxiliaries.h"

using namespace Eigen;

MatrixXd InitialConfig::SetInitialConfiguration() {
	MatrixXd mag(xyz.rows(), 3);
	enum coords : int {x, y, z};
	enum cases : int {
		NotDefined,
		Oblique,
		Homogeneous_x,
		Homogeneous_y,
		Homogeneous_z,
		Homogeneous_mx,
		Homogeneous_my,
		Homogeneous_mz,
		DWall90degXY,
		Random,
		h2h_x,
		h2h_tube,
		Vortex_xy,
		radial
	};

	std::string fileKeyWord = "fromfile_";
	bool readingFromFile = boost::algorithm::starts_with(StartConfiguration, fileKeyWord);
	if (readingFromFile) {
		std::string configFileName = StartConfiguration.erase(0, fileKeyWord.size());
		std::cout << "reading configuration from file: " << configFileName << std::endl;
		VTKReader r(configFileName);
		mag = r.readMag();
		if (!r.checkNodeNumber(xyz.rows())) {
			exit(1);
		}
		startTime = r.getTime();
		return mag;
	}
	StartConfiguration = toLower(StartConfiguration);
	int initialStructure = NotDefined;
	if (StartConfiguration == "oblique") initialStructure = Oblique;
	if (StartConfiguration == "homogeneous_x") initialStructure = Homogeneous_x;
	if (StartConfiguration == "homogeneous_y") initialStructure = Homogeneous_y;
	if (StartConfiguration == "homogeneous_z") initialStructure = Homogeneous_z;
	if (StartConfiguration == "dwall90degxy") initialStructure = DWall90degXY;
	if (StartConfiguration == "homogeneous_mx") initialStructure = Homogeneous_mx;
	if (StartConfiguration == "homogeneous_my") initialStructure = Homogeneous_my;
	if (StartConfiguration == "homogeneous_mz") initialStructure = Homogeneous_mz;
	if (StartConfiguration == "random") initialStructure = Random;
	if (StartConfiguration == "h2h_x") initialStructure = h2h_x;
	if (StartConfiguration == "h2h_tube") initialStructure = h2h_tube;
	if (StartConfiguration == "vortex_xy") initialStructure = Vortex_xy;
	if (StartConfiguration == "radial") initialStructure = radial;

	switch (initialStructure) {
	case Oblique:
		for (unsigned i = 0; i < xyz.rows(); ++i) {
			mag(i, x) = mag(i, y) = mag(i, z) = std::sqrt(1. / 3.);
		}
		break;
	case Homogeneous_x:
		for (unsigned i = 0; i < xyz.rows(); ++i) {
			mag(i, x) = 1.;
			mag(i, y) = mag(i, z) = 0.;
		}
		break;
	case Homogeneous_mx:
			for (unsigned i = 0; i < xyz.rows(); ++i) {
				mag(i, x) = -1.;
				mag(i, y) = mag(i, z) = 0.;
			}
			break;
	case Homogeneous_y:
		for (unsigned i = 0; i < xyz.rows(); ++i) {
			mag(i, y) = 1.;
			mag(i, x) = mag(i, z) = 0.;
		}
		break;
	case Homogeneous_my:
			for (unsigned i = 0; i < xyz.rows(); ++i) {
				mag(i, y) = -1.;
				mag(i, x) = mag(i, z) = 0.;
			}
			break;

	case Homogeneous_z:
	default:
		for (unsigned i = 0; i < xyz.rows(); ++i) {
			mag(i, z) = 1.;
			mag(i, x) = mag(i, y) = 0.;
		}
		break;
	case Homogeneous_mz:
		for (unsigned i = 0; i < xyz.rows(); ++i) {
			mag(i, z) = -1.;
			mag(i, x) = mag(i, y) = 0.;
		}
		break;
	case DWall90degXY: {
		double start = 0.0; // those parameters should be passed to a routine / setter method
		double width = 0.2; // ToDo: replace the hard-coding of these parameters with setters
		double range = 1.;
		double mid = start + range / 2.;
		for (unsigned i = 0; i < xyz.rows(); ++i) {
			double pos = xyz(i, 0);
			double phi = std::tanh((pos - mid) / width) * PhysicalConstants::pi / 4.;
			mag(i, x) = std::cos(phi);
			mag(i, y) = std::sin(phi);
			mag(i, z) = 0.;
		}
	}
		break;
	case h2h_x: {
			int nx = xyz.rows();
			double dw_width = 15.;
			double mid_x = 100.;
			for (int i = 0; i < nx; ++i) {
				mag(i,z) = 0.;
				mag(i,y) = 0.;
				mag(i, x) = 1.;
				double pos_x = xyz(i, x) - mid_x;
				if (pos_x > 0) mag(i, x) = -1.;
				if (std::abs(pos_x) < dw_width / 2. ) {
					mag(i,x) = 0;
					mag(i,y) = 1.;
				}
			}
		}
		break;
	case h2h_tube: {
		int nx = xyz.rows();
		double dw_width = 0.05;
		double mid_z = 1.;
		for (int i = 0; i < nx; ++i) {
			mag(i, x) = 0.;
			mag(i, y) = 0.;
			mag(i, z) = 1.;
			double pos_z = xyz(i, z) - mid_z;
			if (pos_z > 0)
				mag(i, z) = -1.;
			if (std::abs(pos_z) < dw_width / 2.) {
				mag(i, z) = 0;
				mag(i, x) = -xyz(i, y);
				mag(i, y) = xyz(i, x);
				mag.row(i).normalize();
			}
		}
	}
		break;
	case Random: {
		std::random_device rd;
		std::srand(rd());
		int nx = xyz.rows();
		VectorXd theta = PhysicalConstants::pi
				* (VectorXd::Random(nx) + VectorXd::Ones(nx)) / 2.;
		VectorXd phi = PhysicalConstants::pi
				* (VectorXd::Random(nx) + VectorXd::Ones(nx));
		for (int i = 0; i < nx; ++i) {
			mag(i, x) = std::sin(theta(i)) * std::cos(phi(i));
			mag(i, y) = std::sin(theta(i)) * std::sin(phi(i));
			mag(i, z) = std::cos(theta(i));
		}
	}
		break;
	case Vortex_xy: {
		int nx = xyz.rows();
		double core_rad = 2;
		Vector2d center;
		center << 0., 0.;
		Vector2d pos;
		for (int i = 0; i < nx; ++i) {
			pos(0) = xyz(i, x);
			pos(1) = xyz(i, y);
			Vector2d d = pos - center;
			double dist = std::sqrt(d.dot(d));
			if (dist < core_rad) {
				mag(i, x) = 0.;
				mag(i, y) = 0.;
				mag(i, z) = 1.;
			} else {
				mag(i, x) = -d(y);
				mag(i, y) = d(x);
				mag(i, z) = 0;
				mag.row(i).normalize();
			}
		}
	}
		break;
case radial: {
		int nx = xyz.rows();
		double tiny = 1.e-7;
		for (int i = 0; i < nx; ++i) {
			Vector3d pos = xyz.row(i);
			
			if (pos.norm() < tiny) {
				mag(i, x) = 0.;
				mag(i, y) = 0.;
				mag(i, z) = 1.;
			} else {
				mag.row(i) = pos.normalized();
			}
		}
	}
		break;


	}
	return mag;
}


InitialConfig::InitialConfig(const MatrixXd& xyz, std::string StartConfiguration) :
		xyz(xyz), nx(xyz.rows()), StartConfiguration(StartConfiguration), startTime(0.) {
	mag = SetInitialConfiguration();
}


double InitialConfig::GetInitialTime() {
	return startTime;
}


MatrixXd InitialConfig::GetConfiguration() {
	return mag;
}
