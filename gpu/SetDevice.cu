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
 * SetDevice.cu
 *
 *  Created on: Sep 16, 2020
 *      Author: hertel
 */

#include "SetDevice.h"
#include <iostream>

SetDevice::SetDevice(int mydev_) : mydev(mydev_){
	success = false;
	ndev = 0;
}


bool SetDevice::set(){
	success = true;
  cudaError_t err = cudaGetDeviceCount(&ndev);

  if (err != cudaSuccess) {
    std::cout << "No GPU device available. Defaulting to CPU computation" << std::endl;
    success = false;
    return success;
  }
// else {
//    std::cout << "Found " << ndev << " devices." << std::endl;
//  }

  if (mydev > ndev - 1 | mydev < 0) {
    std::cout << "device number " << mydev << " not available. Defaulting to device number 0." << std::endl;
    mydev = 0;
  }
  cudaSetDevice(mydev);
  return success;
}
