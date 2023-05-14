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
 * VTKwriter.h
 *
 *  Created on: Oct 10, 2016
 *      Author: riccardo
 */

#ifndef MESHWRITER_H_
#define MESHWRITER_H_
#include <Eigen/Dense>
#include <string>
#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>
#include <vtkDoubleArray.h>
#include <vtkIntArray.h>
#include <functional>

class MeshWriter {
private:
	Eigen::MatrixXi fel;
	Eigen::MatrixXd xyz;
	std::string name;
	unsigned ntet;
	unsigned nx;
	Eigen::VectorXi NodeMaterial;
	std::string sequenceFileName;
	vtkSmartPointer<vtkUnstructuredGrid> unstructuredGrid;
	void defineUnstructuredVTKGrid();
	std::function <void(MeshWriter*,int, Eigen::MatrixXd&, double, std::string)> writer;
	void VTKwrite(int, Eigen::MatrixXd&, double, std::string);
	void GMVwrite(int, Eigen::MatrixXd&, double, std::string);
	vtkSmartPointer<vtkDoubleArray> setFieldVTK(const Eigen::MatrixXd&, std::string);
	vtkSmartPointer<vtkIntArray> setMaterialVTK(const Eigen::VectorXi& );
	vtkSmartPointer<vtkDoubleArray> setScalarFieldVTK(const Eigen::VectorXd& , std::string );
public:
	void graphicsOutput(int, Eigen::MatrixXd&, double real_time = 0, std::string = "timeInPs");
	void setSequenceFileName(int);
	void selectWriter(std::string);
	void outputVTK(std::string, const Eigen::MatrixXd&, double real_time = 0, std::string = "timeInPs");
	void outputVTK(std::string, const Eigen::MatrixXd&, const Eigen::VectorXd&, std::string = "Magnetization", std::string = "scalar");
	void outputVTK(std::string, const Eigen::MatrixXd&, const Eigen::MatrixXd&, std::string = "Magnetization", std::string = "Field 2");
	void outputGMV(std::string, const Eigen::MatrixXd&, double real_time = 0, std::string = "timeInPs");
	void writeBoundaryGMV(std::string, const Eigen::MatrixXi&, const Eigen::MatrixXd&);
//	void outputVTPolydata(std::string, const Eigen::MatrixXd&, double real_time = 0);
	MeshWriter(){};
	MeshWriter(const Eigen::MatrixXi&, const Eigen::MatrixXd&, std::string, Eigen::VectorXi);
};

#endif /* MESHWRITER_H_ */
