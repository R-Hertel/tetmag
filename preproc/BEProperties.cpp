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

 /* BEProperties.cpp
 *
 *	Created on: Feb. 11, 2018
 *	    Author: riccardo
 *
 * 	Refactored from:
 *
 * TriangularFacet.cpp
 *
 *  Created on: Oct 14, 2016
 *      Author: riccardo
 */

#include "BEProperties.h"
#include <Eigen/Dense>
using namespace Eigen;

TriangleFace::TriangleFace() :
		FacetNodes(0, 0, 0), TetrahedronNumber(0) {}

TriangleFace::TriangleFace(const std::vector<int> FacetNodes_, const int TetrahedronNumber_)
: FacetNodes(Vector3i::Map(FacetNodes_.data())), TetrahedronNumber(TetrahedronNumber_) {}



bool TriangleFace::operator==(const TriangleFace& tf) const {
	if (tf.FacetNodes(0) == this->FacetNodes(0)
			&& tf.FacetNodes(1) == this->FacetNodes(1)
			&& tf.FacetNodes(2) == this->FacetNodes(2))
		return true;
	else
		return false;
}

bool TriangleFace::operator< (const TriangleFace& tf) const {
	if ((this->FacetNodes(0) < tf.FacetNodes(0)) ||
	(this->FacetNodes(0) == tf.FacetNodes(0) && this->FacetNodes(1)  < tf.FacetNodes(1)) ||
	(this->FacetNodes(0) == tf.FacetNodes(0) && this->FacetNodes(1) == tf.FacetNodes(1)	&& this->FacetNodes(2) < tf.FacetNodes(2)))
		return true;
	else
		return false;
}


BEproperties BEproperties::operator=(const TriangleFace& a) {
		this->FacetNodes = a.FacetNodes;
		this->TetrahedronNumber = a.TetrahedronNumber;
		return *this;
}
