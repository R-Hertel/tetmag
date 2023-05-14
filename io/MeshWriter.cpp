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
 * VTKwriter.cpp
 *
 *  Created on: Oct 10, 2016
 *      Author: riccardo
 */
#include "MeshWriter.h"

#include <Eigen/Dense>
#include <string>
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkFieldData.h>
#include <vtkPointData.h>
#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>
#include <vtkUnstructuredGridWriter.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkXMLPolyDataWriter.h>
#include <sstream>
#include <iomanip>
#include <fstream>
using namespace Eigen;

MeshWriter::MeshWriter(const MatrixXi& fel, const MatrixXd& xyz, std::string name, VectorXi Materials) :
		fel(fel), xyz(xyz), name(name) , ntet(fel.rows()), nx(xyz.rows()), NodeMaterial(Materials) {
	defineUnstructuredVTKGrid();
}


void MeshWriter::graphicsOutput(int n,  MatrixXd& vec, double scalar, std::string scalarName) {
	writer(this, n, vec, scalar, scalarName);
}


void MeshWriter::selectWriter(std::string type) {
	if (type == "gmv") {
		writer = &MeshWriter::GMVwrite;
	} else {
		writer = &MeshWriter::VTKwrite;
	}
}


void MeshWriter::defineUnstructuredVTKGrid() {
	unstructuredGrid = vtkSmartPointer <vtkUnstructuredGrid>::New();
	vtkSmartPointer <vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
	points->SetNumberOfPoints(nx);
	for (size_t i = 0; i < nx; ++i) {
		points->SetPoint(i, xyz(i, 0), xyz(i, 1), xyz(i, 2));
	}

	// store finite elements:
	unstructuredGrid->SetPoints(points);
	vtkIdType ptIds[4];
	for (size_t j = 0; j < ntet; j++) {
		for (size_t i = 0; i < 4; i++) {
			ptIds[i] = fel(j, i); // don't try to use Eigen::Map here.
		}
		unstructuredGrid->InsertNextCell(VTK_TETRA, 4, ptIds);
	}
}


void MeshWriter::VTKwrite(int n, MatrixXd& vec, double t, std::string scalarName) {
	setSequenceFileName(n);
	outputVTK(sequenceFileName, vec, t, scalarName);
}


void MeshWriter::setSequenceFileName(int fileNumber) {
	std::stringstream ss;
	ss <<  name << std::setw(7) << std::setfill('0') << fileNumber;
	sequenceFileName = ss.str();
}


void MeshWriter::GMVwrite(int n, MatrixXd& vec, double t, std::string scalarName) {
	setSequenceFileName(n);
	outputGMV(sequenceFileName, vec, t, scalarName);
}


void MeshWriter::outputGMV(std::string filename, const MatrixXd& vec, double t, std::string scalarName ) {
	enum Coords { x, y, z };
	std::ofstream fs;
	fs.open(filename);
	if(!fs) {
		std::cerr << "File " << filename << " could not be opened." << std::endl;
		exit(0);
	}
	fs << "gmvinput ascii\n\n";
	fs << "nodes " << nx << "\n";
	fs << xyz.col(x).transpose() << "\n";
	fs << xyz.col(y).transpose() << "\n";
	fs << xyz.col(z).transpose() << "\n\n";
	fs << "cells " << ntet << "\n\n";
	for (unsigned i = 0; i < ntet; ++i) {
		fs << "tet 4 " << fel.row(i) + RowVector4i::Ones(4) << "\n"; // convert to one-based indexing
	}
	fs << "\n";
	fs << "material 1 1 \nMaterial\n";
	VectorXi mats = NodeMaterial + VectorXi::Ones(nx);
	fs << mats.transpose() << "\n\n";
	fs << "velocity 1\n";
	fs << vec.col(x).transpose() << "\n";
	fs << vec.col(y).transpose() << "\n";
	fs << vec.col(z).transpose() << "\n\n";
//	fs << "probtime \n" << t << std::endl;
	fs << scalarName << "\n" << t << std::endl;
	fs << "endgmv\n";
	fs.close();
}


vtkSmartPointer<vtkDoubleArray> MeshWriter::setFieldVTK(const MatrixXd& vec, std::string fieldName) {
	vtkSmartPointer < vtkDoubleArray > m = vtkSmartPointer < vtkDoubleArray	> ::New();
	m->SetName(fieldName.c_str());
	m->SetNumberOfComponents(3);
	m->SetNumberOfTuples(nx);
	for (size_t i = 0; i < nx; ++i) {
		m->SetTuple3(i, vec(i, 0), vec(i, 1), vec(i, 2));
	}
	return m;
}


vtkSmartPointer<vtkIntArray> MeshWriter::setMaterialVTK(const VectorXi& NodeMaterial) {
	vtkSmartPointer<vtkIntArray> mat = vtkSmartPointer<vtkIntArray>::New();
	mat->SetName("Material");
	mat->SetNumberOfValues(nx);
	for (size_t i = 0; i < nx; ++i) {
		mat->InsertValue(i, NodeMaterial(i));
	}
	return mat;
}


vtkSmartPointer<vtkDoubleArray> MeshWriter::setScalarFieldVTK(const VectorXd& scalarField, std::string fieldName) {
	vtkSmartPointer<vtkDoubleArray> s = vtkSmartPointer<vtkDoubleArray>::New();
	s->SetName(fieldName.c_str());
	s->SetNumberOfValues(nx);
	for (size_t i = 0; i < nx; ++i) {
	  s->InsertValue(i, scalarField(i));
	}
	return s;
}


void MeshWriter::outputVTK(std::string filename, const MatrixXd& vec, double timeValue, std::string scalarName) {
	filename += ".vtu";
	vtkSmartPointer<vtkDoubleArray> m = setFieldVTK(vec, "Magnetization");
	unstructuredGrid->GetPointData()->SetVectors(m);

	// store material
	vtkSmartPointer<vtkIntArray> mat = setMaterialVTK(NodeMaterial);
	unstructuredGrid->GetPointData()->SetScalars(mat);

	// store time data:
//	std::string valueName = "time in ps";
	vtkSmartPointer<vtkDoubleArray> time_ps =  vtkSmartPointer<vtkDoubleArray>::New();
	time_ps->SetNumberOfComponents(1);
	time_ps->SetName(scalarName.c_str());
	time_ps->InsertNextValue(timeValue);
	unstructuredGrid->GetFieldData()->AddArray(time_ps);

	// write file:
	vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
	writer->SetFileName(filename.c_str());
	writer->SetInputData(unstructuredGrid);
	writer->Write();
}


void MeshWriter::outputVTK(std::string filename, const MatrixXd& vec, const VectorXd& scalar, std::string fieldname_1, std::string fieldname_2) {
	filename += ".vtu";

  // store vector field:
  vtkSmartPointer<vtkDoubleArray> m = setFieldVTK(vec, fieldname_1);
  unstructuredGrid->GetPointData()->SetVectors(m);

  // store scalar field
  vtkSmartPointer<vtkDoubleArray> s = setScalarFieldVTK(scalar, fieldname_2);
  unstructuredGrid->GetPointData()->SetScalars(s);

  // write file:
  vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
  writer->SetFileName(filename.c_str());
  writer->SetInputData(unstructuredGrid);
  writer->Write();
}


void MeshWriter::outputVTK(std::string filename, const MatrixXd& vec, const MatrixXd& field, std::string fieldname_1, std::string fieldname_2) { // version with two fields
	filename += ".vtu";
	std::cout << "Writing file: " << filename << std::endl;

	// store magnetization field:
	vtkSmartPointer<vtkDoubleArray> m = setFieldVTK(vec, fieldname_1);
	unstructuredGrid->GetPointData()->SetVectors(m);

	// store material
	vtkSmartPointer < vtkIntArray > mat = setMaterialVTK(NodeMaterial);
	unstructuredGrid->GetPointData()->SetScalars(mat);

	// store additional field
	vtkSmartPointer<vtkDoubleArray> s = setFieldVTK(field, fieldname_2);
	unstructuredGrid->GetPointData()->AddArray(s);

	// write file:
	vtkSmartPointer < vtkXMLUnstructuredGridWriter > writer = vtkSmartPointer < vtkXMLUnstructuredGridWriter > ::New();
	writer->SetFileName(filename.c_str());
	writer->SetInputData(unstructuredGrid);
	writer->Write();
}


void MeshWriter::writeBoundaryGMV(std::string filename, const MatrixXi& bel_b, const MatrixXd& bxyz) {
	enum coords { x, y, z };
	std::ofstream bnd(filename);
	uint bnx = bxyz.rows();
	uint nbel = bel_b.rows();
	if (!bnd) {
		std::cerr << "File " << filename << " could not be opened."	<< std::endl;
		exit(0);
	}
	bnd << "gmvinput ascii\n\n";
	bnd << "nodes " << bnx << "\n";
	bnd << bxyz.col(x).transpose() << "\n";
	bnd << bxyz.col(y).transpose() << "\n";
	bnd << bxyz.col(z).transpose() << "\n\n";
	bnd << "cells " << nbel << "\n\n";
	for (unsigned i = 0; i < nbel; ++i) {
		bnd << "tri 3 " << bel_b.row(i) + RowVector3i::Ones(3) << "\n"; // convert to one-based indexing
	}
	bnd << "\n";
	bnd << "\nendgmv\n";
	bnd.close();
}
