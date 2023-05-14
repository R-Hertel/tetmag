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
 * SetDevice.h
 *
 *  Created on: Sep 16, 2020
 *      Author: hertel
 */

#ifndef GPU_SETDEVICE_H_
#define GPU_SETDEVICE_H_
#include <cuda_runtime_api.h>

class SetDevice {
private:
	int ndev;
	int mydev;
	bool success;
public:
	SetDevice(int);
	bool set();
};

#endif /* GPU_SETDEVICE_H_ */
