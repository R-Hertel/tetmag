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
 * PhysicalConstants.h
 *
 *  Created on: May 29, 2017
 *      Author: riccardo
 */

#ifndef PHYSICALCONSTANTS_H_
#define PHYSICALCONSTANTS_H_
#include <boost/math/constants/constants.hpp>

namespace PhysicalConstants {
	constexpr static double pi = boost::math::constants::pi<double>();
	//double pi = 3.14159265358979323844;
	constexpr static double mu0 = 1.256637061435917295376 * 1.e-6;
	constexpr static double e_charge = 1.6021766208 * 1.e-19;
	constexpr static double e_over_m = -1.758820024 * 1.e11; // C/kg
	constexpr static double g_factor = 2.0023193043618;
	constexpr static double gamma0 = 221276.14886372426114160400; //  rad.m/(A.s)
//	 constexpr static double gamma0 = 1.760859644 * 1.e11; // rad/(s.T) corresponds to the usual 28 [GHz/T] if divided by 2*pi
	constexpr static double mub = 9.274009994 * 1e-24; // Bohr magneton in J/T
}

#endif /* PHYSICALCONSTANTS_H_ */
