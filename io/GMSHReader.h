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
R. Hertel 01.02.2019
 */
#ifndef GMSHREADER_H_
#define GMSHREADER_H_

#include <Eigen/Dense>
#include <string>
#include <vector>

class GMSHReader {
 private:
  Eigen::MatrixXd xyz;
  Eigen::MatrixXi fel;
  std::vector<size_t>  nodeTags;
  size_t nx, nel;
  std::string filename;
  std::vector<bool> isOrphan;
  int hasOrphanNodes();
  void deleteOrphanNodes();
  bool renumberNodes();
  size_t totalOrphans;
 public:
  GMSHReader(std::string);
  Eigen::MatrixXd getNodes();
  Eigen::MatrixXi getElements(); 
  void readPoints();
  void readCells();
  void clean();
  std::vector<size_t> getNodeTags();

};



#endif // GMSHREADER_H_
