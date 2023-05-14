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
 * LindholmCWrapper.cpp
 *
 *  Created on: Mar 5, 2019
 *      Author: riccardo
 */

#include <cmath>
#include <Eigen/Dense>
#include <iostream>
#include "Lindholm.h"
#include "auxiliaries.h"

extern "C"

void LindholmCWrapper (double x0[], double p1[], double p2[], double p3[], double res[]) {
  using namespace Eigen;
  Vector3d Xe(x0);
  Vector3d Pe1(p1);
  Vector3d Pe2(p2);
  Vector3d Pe3(p3);
  double A = triangleArea(Pe1, Pe2, Pe3);
  bool isOrdered = true;
  Vector3d nv = triangleUnitNormVector(Pe1, Pe2, Pe3);
  Lindholm lind;
  Map<Vector3d>(res,3) = lind.weights(Xe, Pe1, Pe2, Pe3, A, nv, isOrdered);
}
