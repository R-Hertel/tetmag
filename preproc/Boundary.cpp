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
 * MeshBoundary.cpp
 *
 *  Created on: Oct 13, 2016
 *      Author: riccardo
 */
#include "Boundary.h"
#include "auxiliaries.h"
#include <algorithm> // for sort
#include <cstddef> // for size_t
#include <boost/math/special_functions/sign.hpp>
#include <queue>
#include <cassert>
#include <iostream>
#include "BEProperties.h"

using namespace Eigen;

MeshBoundary::MeshBoundary(const MatrixXd& xyz_, const MatrixXi& fel_) :
		xyz(xyz_), fel(fel_), nbel(0), bnx(0), totalSurfaceByNodes(0.), totalSurfaceByElements(0.) {}


void MeshBoundary::extractFeatures(){
	findBoundaryElements();
	calculateElementSurface();
	findBoundaryNodes();
	prepareBoundaryNodeIndex();
	calculateBoundaryNodeSurface();
	calculateNormalVector();
	assignNormalVectorsToNodes();
}


void MeshBoundary::prepareBoundaryNodeIndex() {
	uint nx = xyz.rows();
	boundaryNodeIndex = (-1) * VectorXi::Ones(nx); // designed to throw an error if indexing is not correct.
	try {
		for (size_t i = 0; i < bnx; ++i) {
			boundaryNodeIndex[boundaryNodes[i]] = i;
		}
	} catch(...) {
		std::cerr << "indexing error in boundaryNodeIndex" << std::endl;
		exit (EXIT_FAILURE);
	}
}


VectorXd MeshBoundary::getNodeArea() {
	VectorXd ksn = nodeSurfaceArea;
	return ksn;
}


void MeshBoundary::assignNormalVectorsToNodes() {
	uint nx = xyz.rows();
	nv_n = MatrixXd::Zero(nx,3);
	for (uint el = 0; el < nbel; ++el) {
		for (int j = 0; j < 3; ++j) {
			uint currentNode = belprops[el].FacetNodes[j];
			nv_n.row(currentNode) += nv_el.row(el) * belprops[el].TriangleSurface / 3.;
		}
	}
	nv_n.rowwise().normalize();
    nv_n = nv_n.unaryExpr([](double v) { return std::isfinite(v)? v : 0.0; });
}


Vector3d MeshBoundary::CenterOfMass(const MatrixXi& fel, const MatrixXd& xyz, int tet) {
	Vector3d center = Vector3d::Zero();
	for (int node = 0; node < 4; ++node) {
		center += xyz.row(fel(tet, node)) / 4.;
	}
	return center;
}


BEproperties initialize_BE(TriangleFace& f, double& surf) {
	BEproperties BEL;
	BEL = f;
	BEL.TriangleSurface = surf;
	return BEL;
}


void MeshBoundary::calculateNormalVector() {
	nv_el.resize(nbel,3);
	belprops.resize(nbel);
	for (size_t i = 0; i < nbel; ++i) {
		belprops[i] = initialize_BE(BoundaryTriangle[i], ElementSurface[i]);
		belprops[i].TetrahedronCenter = CenterOfMass(fel, xyz, belprops[i].TetrahedronNumber);
		MatrixX3d point(3, 3);
		for (int node = 0; node < 3; ++node) {
			point.row(node) = xyz.row(belprops[i].FacetNodes(node));
		}
		Vector3d n = triangleUnitNormVector(point.row(0), point.row(1), point.row(2));
		Vector3d vectorToCenter = belprops[i].TetrahedronCenter - static_cast<Vector3d>(point.row(0));
		n = -n * boost::math::sign(n.dot(vectorToCenter)); // orient normal vector towards the outside.
		belprops[i].UnitNormalVector = n;
		nv_el.row(i) = n;
		// nv_el.row(i) = n * belprops[i].TriangleSurface; // <-- use this to store normal vectors multiplied by element surface.
	}
/*
// Two simple tests:
   std::cout << "Sums of weighted normal vector components (each should be close to zero):\n";
	double xmean, ymean, zmean;
	xmean = ymean = zmean = 0.;
	for (size_t el = 0; el < nbel; ++el) {
			xmean += nv_el(el, 0) * belprops[el].TriangleSurface ;
			ymean += nv_el(el, 1) * belprops[el].TriangleSurface ;
			zmean += nv_el(el, 2) * belprops[el].TriangleSurface ;
		}
		std::cout << "<x> = " <<  xmean << ", <y> = " << ymean << ", <z> = " << zmean << std::endl;
		xmean = ymean = zmean = 0.;
	std::cout << "Unsigned sum of the normal vectors (total should be entire surface):\n";
	for (size_t el = 0; el < nbel; ++el) {
		xmean += std::abs(nv_el(el, 0)) * belprops[el].TriangleSurface ;
		ymean += std::abs(nv_el(el, 1)) * belprops[el].TriangleSurface ;
		zmean += std::abs(nv_el(el, 2)) * belprops[el].TriangleSurface;
	}
	std::cout << "<x> = " <<  xmean << ", <y> = " << ymean << ", <z> = " << zmean << std::endl;
*/
}


void MeshBoundary::calculateElementSurface() {
	ElementSurface = VectorXd::Zero(nbel);
	for (size_t triangle = 0; triangle < nbel; ++triangle) {
		Vector3i n = BoundaryTriangle[triangle].FacetNodes;
		ElementSurface(triangle) = triangleArea(xyz.row(n(0)), xyz.row(n(1)), xyz.row(n(2)));
	}
	totalSurfaceByElements = ElementSurface.sum();
}


void MeshBoundary::findBoundaryElements() {
	size_t ntet = fel.rows();
	std::vector<TriangleFace> TetrahedronFacets(4 * ntet);
	int FacetNumber = 0;
	std::vector<int> FacetNodes(3);
	for (size_t TetrahedronNumber = 0; TetrahedronNumber < ntet;
			++TetrahedronNumber) {
		for (int j = 0; j < 4; ++j) { //  <-- these nested loops yield the "3 choose 4" combinations.
			int idx = 0;
			for (int k = 0; k < 4; ++k) {
				if (j != k) {
					FacetNodes[idx] = fel(TetrahedronNumber, k);
					++idx;
				}
			}
			std::sort(FacetNodes.begin(), FacetNodes.end());
			TetrahedronFacets[FacetNumber] = TriangleFace(FacetNodes, TetrahedronNumber);
			++FacetNumber;
		}
	}
	std::priority_queue<TriangleFace> pq;
	for (size_t i = 0; i < 4 * ntet; ++i)
		pq.push(TetrahedronFacets[i]); // Filling and automatic sorting based on overloaded operator<

	while (!pq.empty()) {
		TriangleFace BoundaryElementCandidate = pq.top();
		pq.pop();
		if (pq.empty()) { // reached end of queue: triangle is a boundary element.
			BoundaryTriangle.push_back(BoundaryElementCandidate);
			break;
		}
		bool SameAsNextElement = ( BoundaryElementCandidate == pq.top() );
		if (SameAsNextElement)
			pq.pop();
		else
			BoundaryTriangle.push_back(BoundaryElementCandidate);
	}
	nbel = BoundaryTriangle.size();
}


void MeshBoundary::findBoundaryNodes() {
	boundaryNodes.resize(0);
	for (size_t i = 0; i < nbel; ++i)
		for (int j = 0; j < 3; ++j) {
			boundaryNodes.push_back(BoundaryTriangle[i].FacetNodes(j));
		}
	std::sort(boundaryNodes.begin(),boundaryNodes.end()); // <-- stores duplicates consecutively.
	boundaryNodes.erase(unique(boundaryNodes.begin(), boundaryNodes.end()),
			boundaryNodes.end());	// <-- removes consecutive equal node numbers.
	bnx = boundaryNodes.size();

}


void MeshBoundary::calculateBoundaryNodeSurface() {
	int nx = xyz.rows();
	nodeSurfaceArea = VectorXd::Zero(nx); // only non-zero for boundary nodes
	for (size_t CurrentBoundaryElement = 0; CurrentBoundaryElement < nbel; ++CurrentBoundaryElement) {
		for (int j = 0; j < 3; ++j) {
			int CurrentNode = BoundaryTriangle[CurrentBoundaryElement].FacetNodes(j);
			nodeSurfaceArea(CurrentNode) += ElementSurface[CurrentBoundaryElement] / 3.;
		}
	}
	totalSurfaceByNodes = nodeSurfaceArea.sum();
}

/// getter ///
std::vector<TriangleFace> MeshBoundary::getBoundaryElements() {
	return BoundaryTriangle;
}


Eigen::VectorXd MeshBoundary::getElementSurface() {
	return ElementSurface;
}


double MeshBoundary::getTotalSurfaceByElements() {
	return totalSurfaceByElements;
}


std::vector<int> MeshBoundary::getBoundaryNodes() {
	return boundaryNodes;
}


double MeshBoundary::getTotalSurfaceByNodes() {
	return totalSurfaceByNodes;
}


Eigen::MatrixX3d MeshBoundary::getNormalVectors() {
	return nv_el;
}


std::vector<BEproperties> MeshBoundary::getBELproperties() {
	return belprops;
}


Eigen::VectorXi MeshBoundary::getBoundaryNodeIndex() {
	return boundaryNodeIndex;
}


Eigen::MatrixX3d MeshBoundary::getNodeNormalVectors() {
	return nv_n;
}
