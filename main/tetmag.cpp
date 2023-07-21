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

#include <iostream>
#include <Eigen/Dense>
#include <string>
#include "FEMprocessing.h"
#include "SimulationData.h"
#include "TheSimulation.h"
#include "ProgramSpecs.h"
#include "MeshData.h"
#include "Materials.h"
#include <vector>

#include "MeshReaderInterface.h"
#include "MaterialReader.h"
#include "GMSHReader.h"
#include "Timer.h"
#include "auxiliaries.h"
#include "h2interface.h"
#include "tetmagVersion.h"

#ifdef USE_CUDA 
#include "SetDevice.h"
#endif

using namespace Eigen;

int main(int argc, char *argv[]) {
  std::cout << " tetmag -  Copyright (C) 2016-2023 CNRS and Université de Strasbourg\n"
    "Author: Riccardo Hertel\n"
    "This program comes with ABSOLUTELY NO WARRANTY; "
    "This is free software, and you are welcome to redistribute it\n"
    "under certain conditions; type `tetmag --license' for details." << std::endl;
  
  if (argc > 1) {
    std::string argv_s = argv[1];
    if (argv_s == "--license") {
	    printLicense();
	    std::exit(0);
    }
  }
  
	SimulationData sd;
	MeshData msh;
	ProgramSpecs prog;
#ifdef _OPENMP
//#define	EIGEN_DONT_PARALLELIZE
	Eigen::initParallel();
//	Eigen::setNbThreads((omp_get_num_threads() > 2) ? omp_get_num_threads() -1 : 1);
#endif
	printTetmagVersion();

#ifdef __FAST_MATH__ 
std::cerr <<  "-ffast-math is broken, don't use it" << std::endl;
exit(0);
#endif
// read input
	prog.readFile();

	#ifdef USE_CUDA 
		if (prog.solverType == "gpu") {
			SetDevice dev(prog.deviceNumber);
			std::cout << "Using device number " << prog.deviceNumber << std::endl;
			if ( !dev.set() ) {
				prog.solverType = "cpu";
			}
		}
#endif

	if (prog.resuming) {
		std::string lastOutputName = getLastOutputFile(".vtu", prog.getName()+ "0");
		if (lastOutputName != "__NOFILEFOUND__"){
			int digits = 7;
			int lastOutputNumber = extractNumberFromFilename( lastOutputName, digits );
			prog.firstFileNumber = lastOutputNumber + 1;
			prog.initialConfiguration = "fromfile_" + lastOutputName + ".vtu";
			std::cout << "Resuming calculation, starting from output number " << lastOutputNumber << std::endl;
		} else {
			std::cout << "Ignoring resume request: No output files with the name \"" << prog.getName() << "\" found." << std::endl;
			std::cout << "Starting new simulation. " << std::endl;
			prog.resuming = false;
		}
	}
	Timer timer;
	msh.name = prog.getName();
	MaterialReader matReader;
	std::vector<Material> mats = matReader.getMaterials();

//  FEM preprocessing
	MeshReader input(msh.name, prog.meshType);
	input.read();
	msh.xyz = input.getPoints();
	msh.fel = input.getCells();
	PreprocessMesh preproc(msh);
	preproc.prepareAllMatrices();
//	preproc.checkVolumes();
	msh = preproc.getMeshData(); 
	msh.preprocessBEM(prog.useH2);

// problem specification
	msh.ascribeMaterialsToNodes();
	sd.setMaterialParametersAtNodes(msh, mats);
	sd.getProgramData(prog);
	sd.scaleToRealSize();
	msh.selectAnisotropicSurfaces(sd.Ks);
	sd.readLocalFieldProfile();

// simulation
	TheSimulation sim(sd, msh, prog); 
	sim.generateInitialConfiguration();
	timer.start();
	sim.start();
	if (prog.useH2) { deleteHmatrices(); }
	timer.end();
	timer.printDuration();
	return 0;
}

