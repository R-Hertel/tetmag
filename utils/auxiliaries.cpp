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

#include "auxiliaries.h"
#include <limits>
#include <fstream>
#include <string>
#include <iostream>
#include <Eigen/Dense>
#include <iomanip>
#include <cmath>
#include <boost/filesystem.hpp>
#include <sstream>
#include <iostream>
#include <vector>
#include <algorithm>
#include "tetmagVersion.h"

using namespace Eigen;

bool areEqual(double x, double y, int ulp ) { // default value of ulp specified in header file
  double eps = std::numeric_limits<double>::epsilon();
  bool almostEqual = std::abs(x-y) < (eps*std::abs(x+y) * ulp);
  if ((std::abs(x) < eps * ulp) && (std::abs(y) < eps * ulp)) almostEqual = true; // case x,y = 0
  return almostEqual;
}


bool fileExists(std::string fileName) {
	std::ifstream f(fileName.c_str());
	return f.good();
}


void displayProgressBar(double ratio, int barSize) {
	  int barLength = std::round(barSize * ratio);
	  std::cout << "[" << std::left << std::setw(barSize) << std::string(barLength,'#') << "] "
			  << std::round(ratio * 100) << "%\r" ;
}


double triangleArea(const Vector3d& p1, const Vector3d& p2, const Vector3d& p3) {
  Vector3d v12 = p1 - p2;
  Vector3d v13 = p1 - p3;
  double area = 1./2. * (v12.cross(v13)).norm();
  return area;
}


Vector3d triangleUnitNormVector(const Vector3d& p1, const Vector3d& p2, const Vector3d& p3) {
  Vector3d v12 = p1 - p2;
  Vector3d v13 = p1 - p3;
  Vector3d n = (v12.cross(v13)).normalized();
  return n;
}


/* Transfered to inline declaration in header file
 *
MatrixXd cross(const Ref<const MatrixXd>& a, const Ref<const MatrixXd>& b) {
// returns rowwise cross product of two matrices, each with three columns.
	enum coords { x, y, z };
	const int nx = a.rows();
	MatrixXd res(nx, 3);

//	if (a.cols() != 3 || b.cols() != 3 || a.rows() != b.rows()) {
//		std::cout << "Error in auxiliaries::cross(). "
//				"Both matrices must have 3 columns and equal number of rows." << std::endl;
//		exit(1);
//	}

//	assert(a.cols() == 3 && b.cols() == 3 && a.rows() == b.rows());

#ifdef _OPENMP
#pragma omp parallel for shared( a, b, res ) schedule ( runtime )
	for (int i = 0; i < nx; ++i) {
		res(i,x) = a(i,y) * b(i,z) -  a(i,z) * b(i,y);
		res(i,y) = a(i,z) * b(i,x) -  a(i,x) * b(i,z);
		res(i,z) = a(i,x) * b(i,y) -  a(i,y) * b(i,x);
	}
#else
	{
	res.col(x) = a.col(y).cwiseProduct(b.col(z)) - a.col(z).cwiseProduct(b.col(y));
	res.col(y) = a.col(z).cwiseProduct(b.col(x)) - a.col(x).cwiseProduct(b.col(z));
	res.col(z) = a.col(x).cwiseProduct(b.col(y)) - a.col(y).cwiseProduct(b.col(x));
	}
#endif
	return res;
}
*/


std::vector<std::string> getAllOutputFiles(const std::string& ext, std::string basename) {
	namespace fs = ::boost::filesystem;
	std::vector<std::string> theFiles;
	const fs::path& dir = fs::current_path();
	fs::directory_iterator it(dir);
	while (it != fs::end(it)) {
		if (fs::is_regular_file(*it) && it->path().extension() == ext) {
			std::string filename = it->path().filename().string();
			if (filename.find(basename) != std::string::npos) {
				theFiles.push_back(it->path().filename().stem().string());
			}
		}
		++it;
	}
	std::sort(theFiles.begin(), theFiles.end());
	return theFiles;
}


std::string getLastOutputFile(const std::string& ext, std::string basename) {
	std::vector<std::string> theFiles = getAllOutputFiles(ext, basename) ;
	if (theFiles.empty()) {
		return "__NOFILEFOUND__";
	}
	return theFiles.back();
}


int extractNumberFromFilename(const std::string filestem, const int numberOfDigits) {
  std::string numberAsString = filestem.substr(filestem.size() - numberOfDigits);
  std::stringstream numStream(numberAsString);
  int num;
  numStream >> num;
  return num;
}


std::string toLower(std::string s) {
	std::transform(s.begin(), s.end(), s.begin(), [](unsigned char c){return std::tolower(c);});
	return s;
}


void printTetmagVersion() {
	std::cout << "This is tetmag v"	<< TETMAG_VERSION_MAJOR << "." << TETMAG_VERSION_MINOR << "." << TETMAG_VERSION_PATCH << std::endl;
};


void printLicense() {
  std::ifstream f("../License/License.txt");
  std::cout << std::endl;
    if (f.is_open())
        std::cout << f.rdbuf();
}
