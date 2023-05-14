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
 * Hysteresis.h
 *
 *  Created on: Jun 17, 2020
 *      Author: riccardo
 */

#ifndef SPECS_HYSTERESIS_H_
#define SPECS_HYSTERESIS_H_

#include <boost/filesystem.hpp>
namespace fs = boost::filesystem;
#include <fstream>
#include <Eigen/Dense>
#include <string>

class Hysteresis {
private:
	int filenumber;
	double firstField;
	double lastField;
	double deltaH;
	double theta, phi;
	bool increment;
	double fieldNow;
	fs::path mainDirectory;
	fs::path hyteresisOutputDir;	
	void checkSteps();
	void createOutputFolder();
	Eigen::Vector3d fieldDir;
	void setFieldUnitVector();
	void writeHysFileHeader();
public:
	std::ofstream hystLogFile;
	std::ofstream simpleHystFile;
	void init(double, double, double, double, double);
	double nextField();
	bool onBranch();
	void close();	
	void moveFile(std::string);
	std::string getHysFileName(std::string);
	void errorMessage(std::string);
	void writeHysData(const Eigen::VectorXd&, double, const Eigen::MatrixXd& );
};

#endif /* SPECS_HYSTERESIS_H_ */
