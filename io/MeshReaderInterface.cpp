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
 * MeshReader.cpp
 *
 *  Created on: Mar 12, 2019
 *      Author: Riccardo Hertel
 */

#include "MeshReaderInterface.h"

#include <string>
#include <Eigen/Dense>
#include <iostream>
#include "auxiliaries.h"
#include "GMSHReader.h"
#include "VTKReader.h"
#include <numeric> // for std::iota

using namespace Eigen;

MeshReader::MeshReader(std::string name_, std::string type) : fileType(type), name(name_) {

/*
	if (fileType == "msh") {
		fileName = name + ".msh";
	} else if (fileType == "vtk") {
		fileName = name + ".vtk";
	} else if (fileType == "vtu") {
		fileName = name + ".vtu";
*/

	if (fileType == "msh" || fileType == "vtk" || fileType == "vtu") {
		fileName = name + "." + fileType;
	} else {
		std::cerr << "Input file type '" << fileType << "' not supported." << std::endl;
		exit(1);
	}
	if (!fileExists(fileName)) {
		std::cerr << "File " << fileName << " not found." << std::endl;
		exit(0);
	}
}


void MeshReader::read(){
  if (fileType == "msh" ) {  
    GMSHReader theMesh( name );
    theMesh.readPoints();
    theMesh.readCells();
    theMesh.clean();
    xyz = theMesh.getNodes();
    fel = theMesh.getElements();
  } else {
    VTKReader theMesh( fileName );
    theMesh.setFileType( fileType );
    xyz = theMesh.readPoints();
    fel = theMesh.readCells();
  }
}


MatrixXd MeshReader::getPoints() {
  return xyz;
}


MatrixXi MeshReader::getCells() {
  return fel;
}



