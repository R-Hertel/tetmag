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
 * BEMprocessing.cpp
 *
 *  Created on: Oct 11, 2016
 *      Author: riccardo
 */
#include "BEMprocessing.h"
#include "BEProperties.h"
#include "auxiliaries.h"
#include "Boundary.h"
#include "Lindholm.h"
#include "PhysicalConstants.h"
#include "MeshWriter.h"
#include <iostream>
#include <cmath>
#include <vector>
#include <numeric>
#include <string>
#include <ios>
#include <fstream>
#include <map>
#include <iomanip>
#include "h2interface.h"

using namespace Eigen;
double pi = PhysicalConstants::pi;

typedef Matrix<int, Dynamic, Dynamic, RowMajor> rowImat;
typedef Matrix<double, Dynamic, 3, RowMajor> rowDmat;

double getH2err(double sum, int bnx) {
	double errh2 = std::abs((sum - bnx) / bnx);
	return errh2;
}


void displayError(double sum, int bnx) {
	double errh2 = getH2err(sum, bnx);
	double tol = 5.e-3;
	bool ok = ( errh2 < tol );
	if ( ok ) {
//		std::cout <<  "okay." << std::endl; // too much verbosity
	} else {
	std::cout << "H2 Test: Cumulative error = " << std::scientific << errh2  << std::fixed << std::endl;
	std::cerr << "Error: H2 tolerance threshold exceeded." << std::endl;
	exit(1);
	}
}

struct tet_edge {
	tet_edge(int, int);
	int first, second;
	bool isDuplicate;
	int sameAs;
	int number;
	bool operator ==(const tet_edge& e1) {
		if (e1.first == this->first && e1.second == this->second)
			return true;
		else
			return false;
	}
	bool operator <(const tet_edge& b1) const {
		if (b1.first < this->first
				|| (b1.first == this->first && b1.second < this->second))
			return true;
		else
			return false;
	}
	void swap() {
		int temp_p = first;
		first = second;
		second = temp_p;
	}
};


tet_edge::tet_edge(int p1, int p2) :
		first(p1), second(p2), isDuplicate(false), sameAs(-1) {
	if (first > second)
		swap();
}


std::vector<tet_edge> extractRawEdges(int nbel, MatrixXi& bel, MatrixXi& el_ed) {
	std::vector<tet_edge> rawEdges;
	int k = 0;
	for (int i = 0; i < nbel; ++i) {
		for (int j = 0; j < 3; ++j) {
			int p1 = bel(i, j);
			int p2 = bel(i, (j + 1) % 3);
			rawEdges.push_back(tet_edge(p1, p2));
			rawEdges[k].number = k;
			el_ed(i, (j + 2) % 3) = k;
			++k;
		}
	}
	return rawEdges;
}


std::vector<tet_edge> BEMprocessing::prepareEdges() {
	el_ed = MatrixXi::Zero(nbel, 3); // <- will be output
	MatrixXi bel(nbel, 3);
#ifdef _OPENMP
#pragma omp parallel for
#endif
	for (size_t i = 0; i < nbel; ++i) {
		bel.row(i) = belprops[i].FacetNodes; // check global vs surface numbering
	}
	std::vector<tet_edge> rawEdges = extractRawEdges(nbel, bel, el_ed);
//	std::cout << "number of raw edges: " << rawEdges.size() << std::endl;

	std::sort(rawEdges.begin(), rawEdges.end());
	size_t i = 0;

	while (i < rawEdges.size()) {
		if (i + 1 == rawEdges.size())
			break;
		if (rawEdges[i] == rawEdges[i + 1]) {
			rawEdges[i + 1].sameAs = rawEdges[i].number;
			rawEdges[i + 1].isDuplicate = true;
			++i;
		}
		++i;
	}
	std::vector<tet_edge> labeledEdges = rawEdges;
	for (size_t i = 0; i < rawEdges.size(); ++i) {
		labeledEdges[rawEdges[i].number] = rawEdges[i];
	}

	// replace duplicate edges with originals:
#ifdef _OPENMP
#pragma omp parallel for
#endif
	for (size_t i = 0; i < labeledEdges.size(); ++i) {
		if (labeledEdges[i].isDuplicate) {
			for (int j = 0; j < el_ed.size(); ++j) {
				if (el_ed(j) == static_cast<int>(i)) el_ed(j) = labeledEdges[i].sameAs;
			}
		}
	}

//	  std::cout << "\nElement edges after replacing duplicates:\n" << el_ed << std::endl;

	std::map<int, int> imap;
	int newIdx = 0;
	for (size_t i = 0; i < labeledEdges.size(); ++i) {
		if (!labeledEdges[i].isDuplicate) {
			imap.insert(std::make_pair(i, newIdx));
			++newIdx;
		}
	}

	for ( int i = 0; i < el_ed.size(); ++i ) {
		int oldVal = el_ed(i); // <- matrix "unrolled", looping over all elements.
		int newVal = imap.find(oldVal)->second;
		el_ed(i) = newVal;
	}

	std::vector<tet_edge> reindexedEdges; //
	for (size_t i = 0; i < labeledEdges.size(); ++i) {
		if (!labeledEdges[i].isDuplicate) {
			reindexedEdges.push_back(labeledEdges[i]);
		}
	}

#ifndef NDEBUG
	assert(numberOfUniqueEdges == reindexedEdges.size());
	for (int i = 0; i < bel.rows(); ++i) {
		for (int j = 0; j < 3; ++j) {
			int e1 = reindexedEdges[el_ed(i, j)].first;
			int e2 = reindexedEdges[el_ed(i, j)].second;
			assert(bel(i, j) != e1 && bel(i, j) != e2); // edge e[i][j] opposite to node t[i][j]
		}
	}
#endif
//	std::cout << "All edge information is prepared.\n";
	return reindexedEdges;
}


bool isCCW(MatrixXd& xyz, MatrixXi& tri, MatrixXd& nv, int tri_num) {
	Vector3d v1 = xyz.row(tri(tri_num, 0)) - xyz.row(tri(tri_num, 2));
	Vector3d v2 = xyz.row(tri(tri_num, 1)) - xyz.row(tri(tri_num, 0));
	Vector3d nv_temp = v1.cross(v2);
	return (nv_temp.dot(nv.row(tri_num)) > 0);
}

typedef double* pvector;

BEMprocessing::BEMprocessing(const MeshData& msh_) :
		msh(msh_), nx(msh_.xyz.rows()), nel(msh_.fel.rows()), bnx(0), nbel(0) {
}


std::vector<int> BEMprocessing::getBoundaryNodes() {
	return msh.boundaryNodes;
}


MatrixXd BEMprocessing::getLaplaceBEMmatrix() {
	return laplaceBEMmat;
}


void BEMprocessing::reorderNodesCCW() {
	MatrixXi tri(nbel, 3);
	MatrixXd xyz(nx, 3);
	MatrixXd nv(nbel, 3);
	for (size_t i = 0; i < nbel; ++i) {
		tri.row(i) = belprops[i].FacetNodes;
		nv.row(i) =  belprops[i].UnitNormalVector;
	}
	xyz = msh.xyz;
#ifdef _OPENMP	
#pragma	omp parallel for
#endif	
	for (size_t i = 0; i < nbel; ++i) {
#ifdef _OPENMP
# pragma omp critical
#endif		
		{
		bool orderOK = isCCW(xyz, tri, nv, i);
		if (!orderOK) {
			int node_temp = tri(i, 0);
			tri(i, 0) = tri(i, 1);
			tri(i, 1) = node_temp;
			assert(isCCW(xyz, tri, nv, i));
		}
		}
	}
	for (size_t i = 0; i < nbel; ++i) {
		belprops[i].FacetNodes = tri.row(i);
	}
}


void BEMprocessing::H2LibCSuite(std::string name) {
	std::string filename = name + ".h2";
	bool readH2 = fileExists(filename);
	rowDmat nvecs(nbel, 3);
	rowDmat bxy(bnx, 3);
	rowImat tri(nbel, 3);
	VectorXd gram(nbel);
	reorderNodesCCW();
	std::vector<tet_edge> theEdges = prepareEdges();
#ifdef _OPENMP	
#pragma omp parallel for
#endif	
	for (size_t i = 0; i < nbel; ++i) {
#ifdef _OPENMP		
# pragma omp critical
#endif		
		{
			for (int j = 0; j < 3; ++j) {
				tri(i, j) = boundaryNodeIndex[belprops[i].FacetNodes[j]];
			}
			gram(i) = 2. * belprops[i].TriangleSurface;
			nvecs.row(i) = belprops[i].UnitNormalVector;
		}
	}
	for (size_t i = 0; i < bnx; ++i) {
		bxy.row(i) = msh.xyz.row(msh.boundaryNodes[i]);
	}

	rowImat el_ed_R = el_ed;
	int ned = theEdges.size();
	rowImat edgeList(ned, 2);
	for (int i = 0; i < ned; ++i) {
		edgeList(i, 0) = boundaryNodeIndex[theEdges[i].first];
		edgeList(i, 1) = boundaryNodeIndex[theEdges[i].second];
	}
	readH2MatrixFromFile(readH2, filename.c_str());
	defineGeometry(bnx, nbel, ned, bxy.data(), tri.data(), gram.data(),
			nvecs.data(), el_ed_R.data(), edgeList.data());
	pvector result;
	VectorXd u2;
	int bnx2 = getNumberOfVertices();
	VectorXd myOnes = -VectorXd::Ones(bnx2);
	calculateSolidAngles();
	VectorXd diagPart = sub /(4. * pi ) - VectorXd::Ones(sub.size());
	setupH2MatrixC(diagPart.data(), sub.size());
	result = H2_mvp_sub( myOnes.data(), bnx2 );
	u2 = Map<VectorXd>( result, bnx2 );
//	std::cout << "Test of H2 matrix: ";
	displayError(u2.sum(), bnx2);
}


void BEMprocessing::preProcess(bool useH2, std::string name) {

	findBoundaryElementsAndNodes();

	if (!useH2) {
		if (!matrixExists()) {
			calculateSolidAngles();
			assembleLaplaceBEMmatrix();
			storeLaplaceBEMMatrix();
		} else {
			laplaceBEMmat = readLaplaceBEMMatrix();
		}
		testLaplaceMatrix();
	} else {
		H2LibCSuite(name);
	}
//[
//	outputBoundary();
//]
}


void BEMprocessing::outputBoundary() {
	prepareLocalBoundaryMatrices();
	MeshWriter writer;
	writer.writeBoundaryGMV("boundary.gmv", bel_b, bxyz);
}


void BEMprocessing::prepareLocalBoundaryMatrices() {
	bel_b.resize(nbel, 3); // boundary element list with local (bnx) numbering
	bxyz.resize(bnx,3); // coordinates of boundary nodes only.
	for (uint i = 0; i < nbel; ++i) {
		for (int j = 0; j< 3; ++j) {
			int globalNode =  belprops[i].FacetNodes[j];
			bel_b(i,j) = boundaryNodeIndex[globalNode];
		}
	}

	for (uint i = 0; i < bnx; ++i) {
		int globalNode = msh.boundaryNodes[i];
		bxyz.row(i) = msh.xyz.row(globalNode);
	}
}


void BEMprocessing::calculateSolidAngles() {
	sub = VectorXd::Zero(bnx);
	for (size_t i = 0; i < nel; ++i) {
		for (int j = 0; j < 4; ++j) {
			bool isBoundaryNode = (boundaryNodeIndex[msh.fel(i, j)] >= 0);
			if (isBoundaryNode) {
				Vector3d P = msh.xyz.row(msh.fel(i, j));
				Vector3d T1 = msh.xyz.row(msh.fel(i, (j + 1) % 4));
				Vector3d T2 = msh.xyz.row(msh.fel(i, (j + 2) % 4));
				Vector3d T3 = msh.xyz.row(msh.fel(i, (j + 3) % 4));
				sub(boundaryNodeIndex[msh.fel(i, j)]) += lindholm.solidAngle(P, T1, T2, T3);
			}
		}
	}
}


void BEMprocessing::assembleLaplaceBEMmatrix() {
	laplaceBEMmat = MatrixXd::Zero(bnx, bnx);
	VectorXi nds_g(3), nds_b(3);
	std::cout << "Assembling BEM matrix... " << std::endl;
	Matrix3d P;
	for (size_t n = 0; n < bnx; ++n) {
		if (n%20 == 0 ) {
			displayProgressBar( static_cast<double>(n) / static_cast<double>(bnx) );
			std::cout << std::flush;
		}
		int global_n = msh.boundaryNodes[n];
		Vector3d x0 = msh.xyz.row(global_n);
		for (size_t j = 0; j < bnx; ++j) {
			laplaceBEMmat(j, j) = sub(j) / (4. * pi) - 1.;
		}
		for (size_t bel = 0; bel < nbel; ++bel) {
			for (int i = 0; i < 3; ++i) {
				nds_g(i) = belprops[bel].FacetNodes(i);
				P.row(i) = msh.xyz.row(nds_g(i));
			}
			Vector3d nv = belprops[bel].UnitNormalVector;
			double A = belprops[bel].TriangleSurface;
			Vector3d L;
			if ((int)n == boundaryNodeIndex[nds_g(0)] || (int)n == boundaryNodeIndex[nds_g(1)] || (int)n == boundaryNodeIndex[nds_g(2)]) {
				L = Vector3d::Zero(3);
			} else {
				L = lindholm.weights(x0, P.row(0), P.row(1), P.row(2), A, nv);	
			}
			for (int j = 0; j < 3; ++j) {
				int m = boundaryNodeIndex[nds_g(j)];
				laplaceBEMmat(n, m) += L(j);
			}
		}
	}
	displayProgressBar( 100. / 100. );
	std::cout << std::endl;
}


Eigen::MatrixXd BEMprocessing::getBoundaryNodeNormalVectors() {
	return nv_n;
}


void BEMprocessing::findBoundaryElementsAndNodes() {
	MeshBoundary BEmsh(msh.xyz, msh.fel);
	BEmsh.extractFeatures();
	msh.boundaryNodes = BEmsh.getBoundaryNodes();
	belprops = BEmsh.getBELproperties();
	nbel = BEmsh.getBoundaryElements().size();
	bnx = msh.boundaryNodes.size();
	boundaryNodeIndex = BEmsh.getBoundaryNodeIndex();
	nv_n = BEmsh.getNodeNormalVectors();
	nodeArea = BEmsh.getNodeArea(); // only non-zero for boundary nodes.
	std::cout << "Found " << nel << " elements, " << nx << " nodes, and " << bnx << " boundary nodes." << std::endl;
//  std::cout << "Total surface (elements): " << boundary.getTotalSurfaceByElements() << std::endl;
//	std::cout << "Total surface (nodes): " << bem.getTotalSurfaceByNodes() << std::endl;
}


VectorXd BEMprocessing::getNodeArea() {
	return nodeArea;
}


void BEMprocessing::testLaplaceMatrix() {
	std::cout << "Test of BEM matrix: ";
	VectorXd tst = (-1.) * VectorXd::Ones(bnx);
	VectorXd res = laplaceBEMmat * tst;
	displayError(res.sum(), bnx);
}


bool BEMprocessing::matrixExists() {
	bool found = false;
	std::ofstream myFile((msh.name + ".dns").c_str(), std::ios::in);
	if (myFile)
		found = true;
	myFile.close();
	return found;
}


void BEMprocessing::storeLaplaceBEMMatrix() {
	std::string fileName = msh.name + ".dns";
	std::ofstream OutFile(fileName.c_str(), std::ios::out);
	OutFile.write(reinterpret_cast<char*>(&bnx), sizeof(int));
	OutFile.write(reinterpret_cast<char*>(&bnx), sizeof(int));
	OutFile.write(reinterpret_cast<char*>(laplaceBEMmat.data()), sizeof(double) * laplaceBEMmat.size());
	OutFile.close();
}


MatrixXd BEMprocessing::readLaplaceBEMMatrix() {
	std::string fileName = msh.name + ".dns";
	std::ifstream InFile(fileName.c_str(), std::ios::in);
	int nr2 = 0;
	int nc2 = 0;
	InFile.read(reinterpret_cast<char*>(&nr2), sizeof(int));
	InFile.read(reinterpret_cast<char*>(&nc2), sizeof(int));
	assert(nr2 > 0 && nc2 > 0);
	MatrixXd readMat(nr2, nc2);
	InFile.read(reinterpret_cast<char*>(readMat.data()),
			sizeof(double) * nr2 * nc2);
	InFile.close();
	return readMat;
}
