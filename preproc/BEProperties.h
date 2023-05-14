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
    
 /* BEProperties.h
 *
 *	Created on: Feb. 11, 2018
 *	    Author: riccardo
 * 	Refactored from:
 *
 * TriangularFacet.h
 *
 *  Created on: Oct 14, 2016
 *      Author: riccardo
 */

#ifndef BEPROPERTIES_H_
#define BEPROPERTIES_H_
#include <Eigen/Dense>
#include <vector>

struct TriangleFace {
	Eigen::Vector3i FacetNodes;
	int TetrahedronNumber;
	TriangleFace();
	TriangleFace(const std::vector<int>, const int);
	bool operator==(const TriangleFace&) const;
	bool operator< (const TriangleFace&) const;
};

struct BEproperties: public TriangleFace {
	Eigen::Vector3d TetrahedronCenter;
	Eigen::Vector3d UnitNormalVector;
	double TriangleSurface;
	BEproperties operator=(const TriangleFace&); // do I need a "rule of three" destructor and constuctor?
};

#endif /* BEPROPERTIES_H_ */
