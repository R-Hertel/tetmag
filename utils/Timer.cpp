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
 * Timer.cpp
 *
 *  Created on: Jul 10, 2017
 *      Author: riccardo
 */

#include "Timer.h"
#include <iostream>
#include <iomanip>

void Timer::printDuration() {
	std::chrono::hours   hh = std::chrono::duration_cast<std::chrono::hours> (end_tp - start_tp);
	std::chrono::seconds ss = std::chrono::duration_cast<std::chrono::seconds> ((end_tp - start_tp) % std::chrono::minutes(1));
	std::chrono::minutes mm = std::chrono::duration_cast<std::chrono::minutes> ((end_tp - start_tp) % std::chrono::hours(1));
	std::chrono::milliseconds msec = std::chrono::duration_cast<std::chrono::milliseconds> ((end_tp - start_tp) % std::chrono::seconds(1));
	std::cout << "simulation time: " << std::setfill('0')
					<< std::setw(2) << hh.count() << ":"
					<< std::setw(2) << mm.count() << ":"
					<< std::setw(2) << ss.count() << "."
					<< std::setw(3) << msec.count() << std::endl;
}


void Timer::end() {
	end_tp = std::chrono::steady_clock::now();
}


void Timer::start() {
	start_tp = std::chrono::steady_clock::now();
}


double Timer::getRate(double sim_t) {
	double sec = std::chrono::duration<double>(end_tp - start_tp).count();
	return 1000. * sim_t / sec;
}


double Timer::durationInMus() {
  auto mus = std::chrono::duration_cast<std::chrono::microseconds>(total_time);
    return mus.count();
}


void Timer::reset() {
  total_time = std::chrono::microseconds::zero();  
  start();
  end();
}


void Timer::add() {
	end();
	total_time += std::chrono::duration_cast<std::chrono::microseconds>(end_tp - start_tp);
}
