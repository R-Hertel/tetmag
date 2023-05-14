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
 * VTKReader.h
 *
 *  Created on: Oct 6, 2017
 *      Author: riccardo
 */

#ifndef IO_VTKREADER_H_
#define IO_VTKREADER_H_

#include <Eigen/Dense>
#include <string>
#include <functional>
#include <vtkSmartPointer.h>

class VTKReader {
private:
	std::string file;
	std::string fileType;
	double simTime;	
	template <class T> Eigen::MatrixXd readPointsGeneric(vtkSmartPointer<T>);
	template <class T> Eigen::MatrixXi readCellsGeneric(vtkSmartPointer<T>);
	template <class T> Eigen::MatrixXd readMagGeneric(vtkSmartPointer<T>);			
public:
	VTKReader(std::string);
	void setFileType(std::string);
	Eigen::MatrixXd readPoints();
	Eigen::MatrixXi readCells();
	Eigen::MatrixXd readMag();
    Eigen::VectorXi readMaterialType();
	double getTime();
	bool checkNodeNumber(uint);
};

#endif /* IO_VTKREADER_H_ */
