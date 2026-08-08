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
#include <iostream>
#include <fstream>
using namespace Eigen;

MeshWriter::MeshWriter(const MatrixXi& fel, const MatrixXd& xyz, std::string name, VectorXi Materials) :
		fel(fel), xyz(xyz), name(name) , ntet(fel.rows()), nx(xyz.rows()), NodeMaterial(Materials) {
	defineUnstructuredVTKGrid();
}


void MeshWriter::graphicsOutput(int n,  MatrixXd& vec, double scalar, std::string scalarName) {
	VTKwrite(n, vec, scalar, scalarName);
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
	std::string vtuFile = sequenceFileName + ".vtu";
	size_t slashPos = vtuFile.find_last_of('/');
	std::string baseName = (slashPos == std::string::npos) ? vtuFile : vtuFile.substr(slashPos + 1);
	pvdEntries.push_back({t, baseName});
	writePVD();
}


void MeshWriter::writePVD() {
	std::ofstream pvd(name + ".pvd");
	if (!pvd) {
		std::cerr << "Could not open " << name << ".pvd for writing." << std::endl;
		return;
	}
	pvd << "<?xml version=\"1.0\"?>\n";
	pvd << "<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
	pvd << "  <Collection>\n";
	for (const std::pair<double, std::string>& entry : pvdEntries) {
		pvd << "    <DataSet timestep=\"" << std::setprecision(10) << entry.first
		    << "\" part=\"0\" file=\"" << entry.second << "\"/>\n";
	}
	pvd << "  </Collection>\n";
	pvd << "</VTKFile>\n";
}


void MeshWriter::setSequenceFileName(int fileNumber) {
	std::stringstream ss;
	ss <<  name << std::setw(7) << std::setfill('0') << fileNumber;
	sequenceFileName = ss.str();
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

	unstructuredGrid->GetFieldData()->Initialize();
	vtkSmartPointer<vtkDoubleArray> time_ps = vtkSmartPointer<vtkDoubleArray>::New();
	time_ps->SetNumberOfComponents(1);
	time_ps->SetName(scalarName.c_str());
	time_ps->InsertNextValue(timeValue);
	unstructuredGrid->GetFieldData()->AddArray(time_ps);

	// write file:
	vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
	writer->SetFileName(filename.c_str());
	writer->SetInputData(unstructuredGrid);
	writer->SetDataModeToBinary();
	writer->Write();
}


void MeshWriter::addVectorVTK(const MatrixXd& F, std::string fieldName) {
	vtkSmartPointer<vtkDoubleArray> fld = setFieldVTK(F, fieldName);
	unstructuredGrid->GetPointData()->AddArray(fld);
}


void MeshWriter::addScalarVTK(const VectorXd& F, std::string fieldName) {
	vtkSmartPointer<vtkDoubleArray> fld = setScalarFieldVTK(F, fieldName);
	unstructuredGrid->GetPointData()->AddArray(fld);
}


void MeshWriter::closeVTK() {
	std::string filename = name + ".vtu";
	vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
	writer->SetFileName(filename.c_str());
	writer->SetInputData(unstructuredGrid);
	writer->SetDataModeToBinary();
	writer->Write();
	unstructuredGrid->GetPointData()->Initialize();
	unstructuredGrid->GetFieldData()->Initialize();
}


void MeshWriter::writeBoundaryVTK(const MatrixXi& bel_b, const MatrixXd& bxyz) {
	uint bnx = bxyz.rows();
	uint nbel = bel_b.rows();
	vtkSmartPointer<vtkUnstructuredGrid> grid = vtkSmartPointer<vtkUnstructuredGrid>::New();
	vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
	points->SetNumberOfPoints(bnx);
	for (uint i = 0; i < bnx; ++i) {
		points->SetPoint(i, bxyz(i, 0), bxyz(i, 1), bxyz(i, 2));
	}
	grid->SetPoints(points);
	vtkIdType ptIds[3];
	for (uint j = 0; j < nbel; ++j) {
		for (int k = 0; k < 3; ++k) ptIds[k] = bel_b(j, k);
		grid->InsertNextCell(VTK_TRIANGLE, 3, ptIds);
	}
	vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
	writer->SetFileName("boundary.vtu");
	writer->SetInputData(grid);
	writer->SetDataModeToBinary();
	writer->Write();
}


