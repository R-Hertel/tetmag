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
 * MaterialReader.cpp
 *
 *  Created on: May 19, 2017
 *      Author: riccardo
 */

#include "MaterialReader.h"
#include "Materials.h"
#include <sstream>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>
#include "auxiliaries.h"

Material MaterialReader::readMaterialFromFiles(size_t i) {
	std::string file = getMaterialFileName(i+1);
	if (!fileExists(file)) {
		std::cerr << "File " << file << " not found." << std::endl;
		exit(0);
	}
	Material mat;
//	std::cout << "Reading from file " << file << std::endl;
	boost::program_options::variables_map vm;
	namespace opt = boost::program_options;
	opt::options_description desc(
			"material constants provided in materialNNN.dat file:");
	desc.add_options()
	("name",opt::value<std::string>()->default_value(file.substr(0, file.size() - 4)), "Name of material")
	("A", opt::value<double>()->required(), "Exchange constant")
	("Js", opt::value<double>()->required(), "Saturation polarization [T], equal to mu0*Ms")
	("Ku", opt::value<double>()->default_value(0.),	"Uniaxial anisotropy constant")
	("phi_u", opt::value<double>()->default_value(0.), "azimuth angle of uniaxial anisotropy axis")
	("theta_u", opt::value<double>()->default_value(0.), "polar angle of uniaxial anisotropy axis")
	("Kc1",	opt::value<double>()->default_value(0.), "First cubic anisotropy constant")
	("Kc2",	opt::value<double>()->default_value(0.), "Second cubic anisotropy constant")
	("Ks",	opt::value<double>()->default_value(0.), "surface anisotropy constant")
	("phi_Euler", opt::value<double>()->default_value(0.), "Euler angle phi of cubic anisotropy")
	("theta_Euler",	opt::value<double>()->default_value(0.),"Euler angle theta of cubic anisotropy")
	("psi_Euler", opt::value<double>()->default_value(0.),"Euler angle psi of cubic anisotropy")
	("D", opt::value<double>()->default_value(0.),"DMI constant (volume) in J/m2 ");
	std::ifstream cfg_stream(file.c_str());
	try {
		bool allow_unregistered = true; // default boost value is false.
		opt::store(opt::parse_config_file(cfg_stream, desc, allow_unregistered),vm);

	} catch (const opt::error &ex) {
	  std::cout << "ERROR - Input data in material file: ";
	  std::cerr << ex.what() << std::endl;
	  exit(1);
	}
	try {
		notify(vm);
	} catch (boost::program_options::error &e) {
		std::cout << "Error " << e.what() << std::endl;
	}
	if (vm.count("name")) mat.desc =vm["name"].as<std::string>();
	if (vm.count("A")) mat.A = vm["A"].as<double>();
	if (vm.count("Js")) mat.Js = vm["Js"].as<double>();
	if (vm.count("Ku"))	mat.Ku = vm["Ku"].as<double>();
	if (vm.count("phi_u")) mat.phi_u = vm["phi_u"].as<double>();
	if (vm.count("theta_u")) mat.theta_u = vm["theta_u"].as<double>();
	if (vm.count("Kc1")) mat.Kc1 = vm["Kc1"].as<double>();
	if (vm.count("Kc2")) mat.Kc2 = vm["Kc2"].as<double>();
	if (vm.count("Ks")) mat.Ks = vm["Ks"].as<double>();
	if (vm.count("phi_Euler")) mat.phi_Euler = vm["phi_Euler"].as<double>();
	if (vm.count("theta_Euler")) mat.theta_Euler = vm["theta_Euler"].as<double>();
	if (vm.count("psi_Euler")) mat.psi_Euler = vm["psi_Euler"].as<double>();
	if (vm.count("D")) mat.D = vm["D"].as<double>();
	return mat;
}


std::string MaterialReader::getMaterialFileName(size_t num) {
	std::stringstream ss;
	ss << "material" << std::setw(3) << std::setfill('0') << num << ".dat";
	return ss.str();
}


int MaterialReader::countMaterialFiles() {
	std::string lastFile = getLastOutputFile(".dat", "material");
	if (lastFile == "__NOFILEFOUND__"){
		std::cout << "No file containing material data found.\nPlease provide a file named \"material001.dat\"." << std::endl;
		exit(0);
	}
	int num = extractNumberFromFilename(lastFile,3);
	std::cout << "Found " << num << " file" <<
			(num > 1 ? "s" : "") << " with material data." << std::endl;
	return num;
}


MaterialReader::MaterialReader() {
	size_t matTypes = countMaterialFiles();
	if (!matTypes) {
		std::cerr << "Error: Input file of material data not found." << std::endl;
		exit(1);
	}
	for (size_t i = 0; i < matTypes; i++) {
		Material mat = readMaterialFromFiles(i);
		mats.push_back(mat);
	}
	nmats = mats.size();
}


std::vector<Material> MaterialReader::getMaterials() {
	return mats;
}


void MaterialReader::printMaterials() {
	int nmats = mats.size();
	std::cout<< "Total number of materials: " << nmats << "\n";
	for (int i = 0; i < nmats; i++) {
		std::cout << "material #" << i+1 << ":\n";
		std::cout << "====================" << std::endl;
		std::cout << "name = " << mats[i].desc << std::endl;
		std::cout << "A  = " << mats[i].A << std::endl;
		std::cout << "Js = " << mats[i].Js << std::endl;
		std::cout << "Ku = " << mats[i].Ku << std::endl;
		std::cout << "phi_u = " << mats[i].phi_u << std::endl;
		std::cout << "theta_u = " << mats[i].theta_u << std::endl;
		std::cout << "phi_Euler = " << mats[i].phi_Euler << std::endl;
		std::cout << "theta_Euler = " << mats[i].theta_Euler << std::endl;
		std::cout << "psi_Euler = " << mats[i].psi_Euler << std::endl;
		std::cout << "Kc1 = " << mats[i].Kc1 << std::endl;
		std::cout << "Kc2 = " << mats[i].Kc2 << std::endl;
		std::cout << "Ks = " << mats[i].Ks << std::endl;
		std::cout << "D  = " << mats[i].D << std::endl;
		std::cout << "====================" << std::endl;
	}
	exit(0);
}

