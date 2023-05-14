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
 * VTKReader.cpp
 *
 *  Created on: Oct 6, 2017
 *      Author: riccardo
 */

#include <string>
#include <vtkXMLUnstructuredGridReader.h>
#include <vtkUnstructuredGridReader.h>
#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkDoubleArray.h>
#include <vtkCellIterator.h>
#include <iostream>
#include "VTKReader.h"
#include <boost/filesystem/path.hpp>
#include <fstream>
#include "auxiliaries.h"

using namespace Eigen;

typedef vtkSmartPointer <vtkXMLUnstructuredGridReader> vtuPointer;
typedef vtkSmartPointer <vtkUnstructuredGridReader> vtkPointer;


VTKReader::VTKReader(std::string filename_): file(filename_ ){
  boost::filesystem::path p(file);
  std::string ext =  p.extension().string();
  ext.erase(0,1); // remove dot of extension
  simTime = 0;
  setFileType(ext);
}


void VTKReader::setFileType(std::string type) {	
  if (!(type == "vtu" || type == "vtk")) {    
    std::cerr << "VTKReader: file type " << type << " not supported." << std::endl;
    exit(1);
  }
  fileType = type;
  if (!fileExists(file)) {
    std::cerr << "VTKReader:: file " << file << " not found." << std::endl;
    exit(1);
  } 
}


bool VTKReader::checkNodeNumber(uint nx) {
  bool consistency = false;
  vtuPointer reader = vtuPointer::New();
  reader->SetFileName(file.c_str());
  reader->Update();
  unsigned int pointNumber = reader->GetOutput()->GetNumberOfPoints();
  if (pointNumber == nx) { consistency = true; }
  if (!consistency) {
	  std::cerr << "ERROR: While reading file \'" << file << "\'" << std::endl;
	  std::cerr << "Inconsistent number of points. Read " << pointNumber << " but expected " << nx << "." <<std::endl;
  }
  return consistency;
}


MatrixXd VTKReader::readMag() {
  if (fileType == "vtu") {
	  return readMagGeneric(vtuPointer::New());
  } else if (fileType == "vtk") {
	  return readMagGeneric(vtkPointer::New());
  } else {
    std::cerr << "VTKReader: File type " << fileType << " not supported." << std::endl;
    exit(1);
  }
}


MatrixXd VTKReader::readPoints() {
  if (fileType == "vtu") {
	  return readPointsGeneric(vtuPointer::New());
  } else {
	  return readPointsGeneric(vtkPointer::New());
  }
}


MatrixXi VTKReader::readCells() {
  if (fileType == "vtu") {
	return readCellsGeneric(vtuPointer::New());
  } else {
    return readCellsGeneric(vtkPointer::New());
  }
}

template<class T>
MatrixXd VTKReader::readPointsGeneric(vtkSmartPointer<T> reader) {
	reader->SetFileName(file.c_str());
	reader->Update();
	size_t nx = reader->GetOutput()->GetNumberOfPoints();
	MatrixXd RetPoints = MatrixXd::Zero(nx, 3);
	vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints> ::New();
	points->SetNumberOfPoints(nx);
	points = reader->GetOutput()->GetPoints();	
	for (size_t i = 0; i < nx; ++i) {
		for (size_t j = 0; j < 3; j++)	{		
			RetPoints(i, j) = points->GetPoint(i)[j];
		}
	}
	return RetPoints;
}


template<class T>
MatrixXi VTKReader::readCellsGeneric(vtkSmartPointer<T> reader) {
	reader->SetFileName(file.c_str());
	reader->Update();
	int nel = reader->GetOutput()->GetNumberOfCells();
	MatrixXi fel = MatrixXi::Ones(nel, 4);	
	vtkSmartPointer<vtkCellArray> cells = vtkSmartPointer<vtkCellArray> ::New();
	vtkSmartPointer<vtkCellIterator> cellIter = vtkSmartPointer<vtkCellIterator>::Take(reader->GetOutput()->NewCellIterator());
	cellIter->InitTraversal();
	int count = 0;
	for (cellIter->InitTraversal(); !cellIter->IsDoneWithTraversal(); cellIter->GoToNextCell()) {
		for (int j = 0; j < 4; ++j) {
			fel(count, j) = cellIter->GetPointIds()->GetId(j);
		}
		++count;
	}
	return fel;
}


template <class T>
MatrixXd VTKReader::readMagGeneric(vtkSmartPointer<T> reader) {
	reader->SetFileName(file.c_str());
	reader->Update();
	unsigned int nx = reader->GetOutput()->GetNumberOfPoints();
    bool found = reader->GetOutput()->GetPointData()->HasArray("Magnetization");
	if (!found) {
		std::cerr << "No data on the magnetization found in file " << file << std::endl;				
		return MatrixXd::Zero(0,0);
	}
	MatrixXd mag(nx, 3);
	for (uint i = 0; i < nx; ++i) {
       vtkSmartPointer < vtkDataArray > mag_p = reader->GetOutput()->GetPointData()->GetArray("Magnetization");
		for (uint j = 0; j < 3; ++j) {
	       mag(i, j) = mag_p->GetTuple3(i)[j];
		}
	}

	bool fileHasTimeInfo = reader->GetOutput()->GetFieldData()->HasArray("timeInPs");
	if (fileHasTimeInfo) {
		simTime = reader->GetOutput()->GetFieldData()->GetArray("timeInPs")->GetTuple1(0);		
	}

	return mag;
}


double VTKReader::getTime() {
	return simTime;
}


VectorXi VTKReader::readMaterialType() {  
  vtuPointer reader = vtuPointer::New();
  reader->SetFileName(file.c_str());
  reader->Update();
  bool found = reader->GetOutput()->GetPointData()->HasArray("Material");	
  if (!found) {
  	  std::cerr << "No material information found in file " << file << "." << std::endl;
	  return VectorXi::Zero(0);	  
  }
  uint nx = reader->GetOutput()->GetPointData()->GetArray("Material")->GetDataSize();
  VectorXi mats = VectorXi::Zero(nx);
  vtkSmartPointer<vtkDataArray> mat_p = reader->GetOutput()->GetPointData()->GetArray("Material");

  for (uint i = 0; i < nx; ++i) {
	  mats[i] = *mat_p->GetTuple(i);	
  }
  return mats;
}
