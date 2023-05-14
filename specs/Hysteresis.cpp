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
 * Hysteresis.cpp
 *
 *  Created on: Jun 17, 2020
 *      Author: riccardo
 */

#include "Hysteresis.h"
#include <sstream>
#include <iostream>
#include <cmath>
#include <cstdlib>
#include <string>
#include <iomanip>

#include "PhysicalConstants.h"


void Hysteresis::init(const double firstField_, const double lastField_, double deltaH_, double theta_, double phi_) {
	filenumber = 0;
	firstField = firstField_;
	lastField = lastField_;
	deltaH = deltaH_;
	theta = theta_;
	phi = phi_;
	increment = true;
	fieldNow = firstField;
	setFieldUnitVector();
	checkSteps();
	createOutputFolder();
	writeHysFileHeader();
}


void Hysteresis::writeHysFileHeader() {
	simpleHystFile << "# Data in Columns: (1) Hysteresis field value (magnitude) in T "
			<< "(2) Projection of <m> along the applied field direction." << std::endl;
}


void Hysteresis::setFieldUnitVector() {
	enum {x, y, z};
	double deg2rad = PhysicalConstants::pi / 180.;
	fieldDir[x] = std::sin(theta * deg2rad) * std::cos(phi * deg2rad);
	fieldDir[y] = std::sin(theta * deg2rad) * std::sin(phi * deg2rad);
	fieldDir[z] = std::cos(theta * deg2rad) ;
}


void Hysteresis::checkSteps() {
  deltaH = std::abs(deltaH); // sign of increment/decrement depends on first and last field value
  double tol = 1.e-8;
  if ( std::abs( firstField - lastField ) < tol ) {
    std::cerr << "ERROR (HYS): Difference between initial and final field value is too small." << std::endl;
    exit(1);
  }
  if ( deltaH > std::abs( firstField - lastField ) ) {
    std::cerr << "ERROR (HYS): Field step larger than difference between first and last field value." << std::endl;
    exit(1);
  }
  if ( deltaH == 0 ) {
    std::cerr  << "ERROR (HYS): Increment deltaH is zero: "  << deltaH << ".\n";
    exit(1);
  }
  if ( firstField > lastField ) {
    increment = false;
    deltaH = -deltaH;
  }
}


double Hysteresis::nextField(){
	fieldNow += deltaH;
	return fieldNow;
}


bool Hysteresis::onBranch() {
	 return (( increment && ( fieldNow <= lastField )) || ( !increment && ( fieldNow >= lastField )));
}


void Hysteresis::createOutputFolder() {
	mainDirectory = fs::current_path();
	std::string dstFolder = "Hysteresis_";
	std::string pathName(dstFolder);

	for (int i = 1; fs::exists(pathName) && i < 10; ++i) {
		std::stringstream ss;
		ss << dstFolder << i << "_";
		pathName = ss.str();
	}

	if (fs::exists(pathName)) {
		std::cerr << "Could not create folder for hysteresis data."
				<< std::endl;
		exit(1);
	}

	boost::filesystem::create_directory(pathName);
	fs::current_path(fs::system_complete(pathName.c_str()));
	fs::path new_full_path(fs::current_path());
	hyteresisOutputDir = new_full_path;
	hystLogFile.open("hysteresis.log");
	simpleHystFile.open("simpleHyst.dat");
	fs::current_path(mainDirectory);
}


void Hysteresis::moveFile(std::string fileName) {
	fs::current_path(mainDirectory);
	fs::path src(fileName);
	fs::rename(src, hyteresisOutputDir / fileName);
}


std::string Hysteresis::getHysFileName(std::string name) {
	std::stringstream ss;
	ss << name << std::setw(3) << std::setfill('0') << filenumber;
	std::string fileName = ss.str();
	filenumber++;
	return fileName;
}


void Hysteresis::errorMessage(std::string type) {
	std::cerr
			<< "ERROR: Time-varying features are not allowed in hysteresis loop calculations.\n"
			<< "Please unset the option '" << type << "' and try again."
			<< std::endl;
	exit(1);
}


void Hysteresis::writeHysData(const Eigen::VectorXd& nodeVolume, double totalVolume, const Eigen::MatrixXd& mag) {
	Eigen::VectorXd magAlongH = (mag * fieldDir).rowwise().sum();
	double meanMtimesH = nodeVolume.dot(magAlongH) / totalVolume;
	simpleHystFile << fieldNow / 1000. << "\t" << meanMtimesH << std::endl;
}


void Hysteresis::close() {
	hystLogFile.close();
	simpleHystFile.close();
}
