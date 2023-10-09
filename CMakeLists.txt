#    tetmag - A general-purpose finite-element micromagnetic simulation software package
#    Copyright (C) 2016-2023 CNRS and Université de Strasbourg
#
#    Author: Riccardo Hertel
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU Affero General Public License as
#    published by the Free Software Foundation, either version 3 of the
#    License, or (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU Affero General Public License for more details. 
# 
#    Contact: Riccardo Hertel, IPCMS Strasbourg, 23 rue du Loess, 
#    	     67034 Strasbourg, France.
#	     riccardo.hertel@ipcms.unistra.fr
#	     
#    You should have received a copy of the GNU Affero General Public License
#    along with this program.  If not, see <https://www.gnu.org/licenses/>.
#
cmake_minimum_required(VERSION 3.13)
set(CMAKE_CXX_STANDARD 14)

project(tetmag VERSION 1.7.6 LANGUAGES C CXX)
# These options can be set with ccmake or with cmake -D 
option(TETMAG_USE_OPENMP "Compile with OpenMP" ON)
option(TETMAG_USE_CUDA "Enable CUDA" ON)
option(TETMAG_FIND_CUDA_ARCH "Autodetect available GPU architectures" ON)
option(TETMAG_USE_CVODE_5 "Use CVODE v5 instead of latest release" OFF)
set(TETMAG_SET_CUDA_ARCHS "61;75" CACHE STRING "manually set GPU architectures")
set(TETMAG_HOST_COMPILER "/usr/bin/g++" CACHE FILEPATH "CUDA-compatible host compiler") 
# ^^^ select host compiler version that is compatible with the installed CUDA version.
# (see https://gist.github.com/ax3l/9489132 )


include(CheckLanguage)
check_language(CUDA)
if (CMAKE_CUDA_COMPILER)
  set(CMAKE_CUDA_HOST_COMPILER ${TETMAG_HOST_COMPILER} )  
  enable_language(CUDA)
endif(CMAKE_CUDA_COMPILER)

if(TETMAG_USE_CUDA)  
  include(FindCUDA) 
  if (CUDA_FOUND)
    set(USE_CUDA ON)
    if (CMAKE_CUDA_COMPILER_VERSION LESS 10.1)    
      message("-- A CUDA version of at least 10.1 is needed (found unsuitable version ${CMAKE_CUDA_COMPILER_VERSION})")
      message("-- No suitable CUDA version found. Not using GPU.")
      set(USE_CUDA OFF)     
    endif()
  else()
    message("-- CUDA not found. Not using GPU.")
    set(USE_CUDA OFF)
  endif(CUDA_FOUND)
else()
  set(USE_CUDA OFF)  
endif(TETMAG_USE_CUDA)

if (TETMAG_USE_OPENMP)
  include(FindOpenMP)  
  if (NOT OpenMP_FOUND)
    message("-- OpenMP not found.")
  endif(NOT OpenMP_FOUND)
  set(USE_OpenMP ${OpenMP_FOUND})
else()
  set(USE_OpenMP OFF)
endif(TETMAG_USE_OPENMP)

message("-- Build with options CUDA=${USE_CUDA} and OpenMP=${USE_OpenMP}")

set(EXTDIR ${PROJECT_SOURCE_DIR}/external )
set(LIBDIR ${PROJECT_SOURCE_DIR}/lib )

if (APPLE)
  link_directories(/usr/local/lib)
  link_directories(${LIBDIR})
  include_directories(/usr/local/include)
endif()

#set_property(GLOBAL PROPERTY RULE_MESSAGES OFF)
#set(CMAKE_VERBOSE_MAKEFILE OFF) 


set(default_build_type "Release")

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "${default_build_type}")
endif()

message("-- Build type set to '${CMAKE_BUILD_TYPE}'.")

configure_file(${PROJECT_SOURCE_DIR}/main/tetmagVersion.h.in ${PROJECT_SOURCE_DIR}/main/tetmagVersion.h)

#string(APPEND CMAKE_CXX_FLAGS_RELEASE " -DEIGEN_USE_BLAS")
#string(APPEND CMAKE_CXX_FLAGS_RELEASE " -pg")

if (TETMAG_USE_CVODE_5)
  add_compile_definitions(USE_CVODE_5)
endif() 

if (CMAKE_BUILD_TYPE MATCHES Release)
  #string(APPEND CMAKE_CXX_FLAGS_RELEASE " -mavx " ) # not compatible with CUDA
  #string(APPEND CMAKE_CXX_FLAGS_RELEASE " -march=native") # not compatible with CUDA
  string(APPEND CMAKE_CXX_FLAGS_RELEASE " -fno-math-errno -Wno-unused -Wall -DNDEBUG" )
  message("-- Build flags: ${CMAKE_CXX_FLAGS_RELEASE} ")  
endif()

if (CMAKE_BUILD_TYPE MATCHES Debug)
	string(APPEND CMAKE_CXX_FLAGS_DEBUG " -ggdb3")
	message("-- Build flags: ${CMAKE_CXX_FLAGS_DEBUG} ")
endif()

find_package(Eigen3 3.3.7 REQUIRED NO_MODULE)
find_package(Boost REQUIRED COMPONENTS filesystem program_options system thread)
find_package(Threads REQUIRED)
find_package(VTK REQUIRED)
find_library(NCDF netcdf REQUIRED)
add_definitions(-DUSE_NETCDF)

find_package(BLAS)
find_package(LAPACK)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/lib/)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PROJECT_SOURCE_DIR}/cmake/)

set(PREPROC ${PROJECT_SOURCE_DIR}/preproc/)
set(UTILS ${PROJECT_SOURCE_DIR}/utils/)
set(IO ${PROJECT_SOURCE_DIR}/io/)
set(SPECS ${PROJECT_SOURCE_DIR}/specs/)
set(MAIN ${PROJECT_SOURCE_DIR}/main/)
set(GPU ${PROJECT_SOURCE_DIR}/gpu/)

if (USE_CUDA)   	
  include(${GPU}/CMakeLists.txt)
endif(USE_CUDA)

include(${EXTDIR}/CMakeLists.txt)   # install external dependencies

include_directories(${PREPROC} ${UTILS} ${IO} ${SPECS} ${MAIN} ${GPU} ${EIGEN3_INCLUDE_DIR} ${VTK_INCLUDE_DIRS})

if (VTK_VERSION VERSION_LESS "8.90.0")
  include(${VTK_USE_FILE})
endif() 

set (PREPROC_FILES
  ${PREPROC}BEMprocessing.cpp  
  ${PREPROC}BEProperties.cpp
  ${PREPROC}Boundary.cpp
  ${PREPROC}FEMprocessing.cpp
  ${PREPROC}Lindholm.cpp
  ${PREPROC}LindholmCWrapper.cpp  
)

file(GLOB PREPROC_SRC ${PREPROC_FILES})
add_library(prep STATIC ${PREPROC_SRC})
set(PREPLIB prep)
target_link_libraries(${PREPLIB} ${VTK_LIBRARIES})

if (NOT VTK_VERSION VERSION_LESS "8.90.0")  
  vtk_module_autoinit(TARGETS ${PREPLIB} MODULES ${VTK_LIBRARIES})
endif()

set (TETMAG_FILES
  ${MAIN}tetmag.cpp  
  ${MAIN}AMGCLSolver.cpp
  ${UTILS}auxiliaries.cpp
  ${MAIN}DemagField.cpp
  ${MAIN}EffFieldCalc.cpp
  ${IO}GMSHReader.cpp
  ${MAIN}h2interface.c
  ${SPECS}Hysteresis.cpp
  ${SPECS}InitMag.cpp    
  ${IO}MaterialReader.cpp
  ${SPECS}Materials.cpp
  ${SPECS}MeshData.cpp 
  ${IO}MeshReaderInterface.cpp
  ${IO}MeshWriter.cpp
  ${UTILS}PhysicalConstants.cpp
  ${SPECS}ProgramSpecs.cpp
  ${UTILS}ProcessInformation.cpp
  ${SPECS}SimulationData.cpp
  ${MAIN}SolverFactory.cpp
  ${MAIN}TheLLG.cpp
  ${MAIN}TheSimulation.cpp
  ${UTILS}Timer.cpp
  ${IO}VTKReader.cpp
)
 
if (USE_CUDA)
  set (CVODE_CUDA sundials_nveccuda)
	cuda_add_executable(${PROJECT_NAME} ${TETMAG_FILES} ${CUDA_FILES})
	target_link_libraries(${PROJECT_NAME} ${CUDA_cusparse_LIBRARY} ${CVODE_CUDA})	
else()
 	 add_executable(${PROJECT_NAME} ${TETMAG_FILES} )
endif(USE_CUDA)	

add_dependencies(${PROJECT_NAME} amgcl sundials_cvode gmsh h2lib)
target_link_libraries(${PROJECT_NAME} Eigen3::Eigen)
target_link_libraries(${PROJECT_NAME} Threads::Threads)
target_link_libraries(${PROJECT_NAME} Boost::boost Boost::filesystem Boost::program_options Boost::system Boost::thread)
target_link_libraries(${PROJECT_NAME} ${H2LIB} ${GMSHLIB} ${CVODELIB})
target_link_libraries(${PROJECT_NAME} ${VTK_LIBRARIES})
target_link_libraries(${PROJECT_NAME} ${PREPLIB} ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES} ${NCDF})

if(USE_OpenMP)
 	target_link_libraries( ${PROJECT_NAME} OpenMP::OpenMP_CXX ) 
  set (CVODE_OMP sundials_nvecopenmp)
	add_dependencies(${PROJECT_NAME} ${CVODE_OMP})
	target_link_libraries(${PROJECT_NAME} ${CVODE_OMP})
endif(USE_OpenMP)	

set_target_properties(${PROJECT_NAME} PROPERTIES LINKER_LANGUAGE CXX)

target_include_directories(${PROJECT_NAME} PUBLIC ${EIGEN3_INCLUDE_DIR})
target_include_directories(${PROJECT_NAME} PUBLIC $<BUILD_INTERFACE:${EXTDIR}/gmsh/include>)
target_include_directories(${PROJECT_NAME} PUBLIC $<BUILD_INTERFACE:${EXTDIR}/h2lib/Library>)
target_include_directories(${PROJECT_NAME} PUBLIC $<BUILD_INTERFACE:${amgcl_SOURCE_DIR}>)

if (NOT VTK_VERSION VERSION_LESS "8.90.0")  
  vtk_module_autoinit(TARGETS ${PROJECT_NAME} MODULES ${VTK_LIBRARIES})
endif()
