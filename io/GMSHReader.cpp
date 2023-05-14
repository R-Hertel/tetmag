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
 * GMSHReader.cpp
 *
 *  Created on: Feb 1, 2019
 *      Author: Riccardo Hertel 
 *
 */
#include <Eigen/Dense>
#include <string>
#include <gmsh.h>
#include <vector>
#include <map>
#include <iostream>
#include "GMSHReader.h"
#include <numeric> // for std::accumulate
using namespace Eigen;
  
GMSHReader::GMSHReader(std::string name) : filename(name + ".msh") {
  gmsh::initialize();
  gmsh::open(filename);  
}


void GMSHReader::readPoints() {
  std::vector<double> coords;
  std::vector<double> parametricCoord; // not used here. Only required to match the argument list of getNodes()
  int dim = -1; // negative value: get ALL nodes, not only those belonging to elements of specific dim
  int tag = -1; // negative value: get get all nodes, irrespective of tag
  const bool includeBoundary = true;	   
  gmsh::model::mesh::getNodes( nodeTags, coords, parametricCoord, dim, tag, includeBoundary ); 
  nx =  nodeTags.size();
  xyz.resize(nx,3);
  xyz = Map<Matrix<double, Dynamic, Dynamic, RowMajor> >(coords.data(), nx, 3);
  nodeTags = getNodeTags();
}


void GMSHReader::readCells() {
  int dim = 3; // get only three-dimensional elements (tetrahedrons)
  int tag = -1; // negative: get get all nodes, irrespective of tag
  std::vector<int> elementTypes; // unnecessary if mesh contains only one type of elements. Needed as argument of getElements()
  std::vector<std::vector<size_t> >  connectivities;
  std::vector<std::vector<size_t> >  elementTags;
  gmsh::model::mesh::getElements( elementTypes, elementTags, connectivities, dim, tag );
  nel =  elementTags[0].size();
  fel.resize(nel,4);
  std::vector<int> myc(nel *  4);
  for (size_t i = 0; i < nel*4; ++i)
      myc[i] = static_cast<int>(connectivities[0][i]); // Map currently does not work with size_t type
  fel = Map<Matrix<int, Dynamic, Dynamic, RowMajor> >(myc.data(), nel, 4);
}


bool GMSHReader::renumberNodes(){
  bool hasFailed = false;
  std::map<int,int> nodeTagsMap;
  for (size_t i = 0; i < nx; ++i) {
    nodeTagsMap.insert(std::make_pair(nodeTags[i],i));
  }
  for (size_t i = 0; i < nel; ++i) {
    for (size_t j = 0; j < 4; ++j) {
      std::map<int,int>::iterator posKey = nodeTagsMap.find(fel(i,j));
      if (posKey != nodeTagsMap.end()) {
	fel(i,j) = posKey->second;
      } else { hasFailed = true; }
    }
  }
  return hasFailed;
}


void GMSHReader::clean() {
  renumberNodes(); 
  if ( hasOrphanNodes() ) {
    deleteOrphanNodes();
    renumberNodes();
  }
};


int GMSHReader::hasOrphanNodes() {
  isOrphan.resize( nx );
  isOrphan.assign( nx, true );
  for (size_t i = 0; i < nel; ++i) {
    for (int j = 0; j < 4; ++j) {
      isOrphan[fel(i,j)] = false;
    }
  }
  totalOrphans = std::accumulate(isOrphan.begin(), isOrphan.end(), 0);
  return totalOrphans;
}


void GMSHReader::deleteOrphanNodes() {
  MatrixXd xyz_new(nx - totalOrphans, 3);
  std::cout << "Removing " << totalOrphans << " orphan node" << (totalOrphans >1 ? "s" : "") << std::endl;
  int reindex = 0;
  nodeTags.resize(nx - totalOrphans);
  for (size_t i = 0; i < nx; ++i) {
    if (!isOrphan[i]) {
      xyz_new.row(reindex) = xyz.row(i);     
      nodeTags[reindex] = i;
      ++reindex;
    }
  }
  xyz = xyz_new;
};


Eigen::MatrixXi GMSHReader::getElements() {
  return fel;
}


Eigen::MatrixXd GMSHReader::getNodes() {
  return xyz;
}


std::vector<size_t> GMSHReader::getNodeTags() {
  return nodeTags;
}

