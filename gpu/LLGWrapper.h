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
 * LLGWrapper.h
 *
 *  Created on: Sep 23, 2019
 *      Author: riccardo
 */

#ifndef GPU_LLGWRAPPER_H_ 
#define GPU_LLGWRAPPER_H_ 

#include <thrust/device_vector.h>
//#include <boost/numeric/odeint.hpp>
//#include <boost/numeric/odeint/external/thrust/thrust.hpp>
#include "TheLLG.h"
#include <vector>
#include <memory>

////  sundials
#include <sunlinsol/sunlinsol_spgmr.h>  //access to SPGMR SUNLinearSolver
#include <cvode/cvode_spils.h> // access to CVSpils interface
#include <cvode/cvode.h> // prototypes for CVODE fcts., consts.
#include <nvector/nvector_serial.h>  // access to serial N_Vector
#include <sundials/sundials_types.h>  // defs. of realtype, sunindextype
#include <sunlinsol/sunlinsol_dense.h> /* access to dense SUNLinearSolver          */
#include <sunmatrix/sunmatrix_dense.h> /* access to dense SUNMatrix                */
#include <nvector/nvector_cuda.h>
#include <thrust/device_ptr.h>

#include <memory>

typedef thrust::device_vector< double > dev_vec;

class LLGWrapper {
	int nx3;
public:
//	void setMag(std::vector<double>&);
	void init(int);
	~LLGWrapper();
};


#endif /* GPU_LLGWRAPPER_H_ */
