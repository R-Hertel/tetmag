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
 * MaterialReader.h
 *
 *  Created on: May 19, 2017
 *      Author: riccardo
 */

#ifndef MATERIALREADER_H_
#define MATERIALREADER_H_
#include "Materials.h"
#include <vector>
#include <string>

class MaterialReader {
private:
	std::string name;
	std::vector<Material> mats;
	int nmats;
	std::string getMaterialFileName(size_t);
	Material readMaterialFromFiles(size_t);
	int countMaterialFiles();
public:
	std::vector<Material> getMaterials();
	MaterialReader();
	void printMaterials();
};

#endif /* MATERIALREADER_H_ */
