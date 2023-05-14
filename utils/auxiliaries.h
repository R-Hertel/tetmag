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
 * auxiliaries.hpp
 *
 *  Created on: Feb 23, 2017
 *      Author: riccardo
 */

#ifndef AUXILIARIES_H_
#define AUXILIARIES_H_

#include <Eigen/Dense>
#include <string>
#include <vector>

bool areEqual(double , double , int ulp = 100) ; // ulp: units in the last place; tolerance.

Eigen::Vector3d triangleUnitNormVector(const Eigen::Vector3d& p1, const Eigen::Vector3d& p2, const Eigen::Vector3d& p3);

double triangleArea(const Eigen::Vector3d&, const Eigen::Vector3d&, const Eigen::Vector3d& );

bool fileExists(std::string fileName);

void displayProgressBar(double, int barSize = 40);

inline Eigen::MatrixXd cross(const Eigen::Ref<const Eigen::MatrixXd>& a, const Eigen::Ref<const Eigen::MatrixXd>& b) {
// returns rowwise cross product of two matrices, each with three columns.
	enum coords { x, y, z };
	const int nx = a.rows();
	Eigen::MatrixXd res(nx, 3);
//	assert(a.cols() == 3 && b.cols() == 3 && a.rows() == b.rows());


#ifdef _OPENMP
#pragma omp parallel for schedule ( runtime )
//#pragma omp parallel for shared( res ) schedule ( runtime )
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


std::string getLastOutputFile(const std::string& , std::string );

std::vector<std::string> getAllOutputFiles(const std::string& , std::string );

int extractNumberFromFilename(const std::string, const int);

std::string toLower(std::string);

void printTetmagVersion();

void printLicense();

#endif /* AUXILIARIES_H_ */

