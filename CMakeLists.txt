#    tetmag - A general-purpose finite-element micromagnetic simulation software package
#    Copyright (C) 2016-2026 CNRS and Université de Strasbourg
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
cmake_minimum_required(VERSION 3.18)

project(tetmag VERSION 1.9.2 LANGUAGES C CXX)
# These options can be set with ccmake or with cmake -D 
option(TETMAG_USE_OPENMP "Compile with OpenMP" ON)
option(TETMAG_USE_CUDA "Enable CUDA" ON)
option(TETMAG_FIND_CUDA_ARCH "Autodetect available GPU architectures" ON)
set(TETMAG_SET_CUDA_ARCHS "61;75" CACHE STRING "manually set GPU architectures")
set(TETMAG_HOST_COMPILER "/usr/bin/g++" CACHE FILEPATH "CUDA-compatible host compiler")
# ^^^ select host compiler version compatible with the installed CUDA version.
# (see https://gist.github.com/ax3l/9489132 )

# Installation prefix
if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  set(CMAKE_INSTALL_PREFIX "/usr/local" CACHE PATH "Installation prefix" FORCE)
endif()
set(INSTALL_BIN_DIR "bin" CACHE PATH "")
set(INSTALL_LIB_SUBDIR "lib/tetmag" CACHE PATH "")


set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include(CheckLanguage)
check_language(CUDA)
if(CMAKE_CUDA_COMPILER)
  set(CMAKE_CUDA_STANDARD 17)
  set(CMAKE_CUDA_HOST_COMPILER ${TETMAG_HOST_COMPILER})
  enable_language(CUDA)
endif(CMAKE_CUDA_COMPILER)

if(TETMAG_USE_CUDA)
  find_package(CUDAToolkit QUIET)
  if(CUDAToolkit_FOUND)
    set(USE_CUDA ON)
    if(CMAKE_CUDA_COMPILER_VERSION LESS 10.1)
      message(STATUS "A CUDA version of at least 10.1 is needed (found unsuitable version ${CMAKE_CUDA_COMPILER_VERSION})")
      message(WARNING "No suitable CUDA version found. Not using GPU.")
      set(USE_CUDA OFF)
    endif()
  else()
    message(WARNING "CUDA not found. Not using GPU.")
    set(USE_CUDA OFF)
  endif(CUDAToolkit_FOUND)
else()
  set(USE_CUDA OFF)
endif(TETMAG_USE_CUDA)

if(TETMAG_USE_OPENMP)
  include(FindOpenMP)
  if(NOT OpenMP_FOUND)
    message(STATUS "OpenMP not found.")
  endif(NOT OpenMP_FOUND)
  set(USE_OpenMP ${OpenMP_FOUND})
else()
  set(USE_OpenMP OFF)
endif(TETMAG_USE_OPENMP)

message(STATUS "Build with options CUDA=${USE_CUDA} and OpenMP=${USE_OpenMP}")

set(EXTDIR ${PROJECT_SOURCE_DIR}/external)
set(LIBDIR ${PROJECT_SOURCE_DIR}/lib)

if(APPLE)
  link_directories(/usr/local/lib)
  link_directories(${LIBDIR})
  include_directories(/usr/local/include)
endif()

set(default_build_type "Release")

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "${default_build_type}")
endif()

message(STATUS "Build type set to '${CMAKE_BUILD_TYPE}'.")

configure_file(${PROJECT_SOURCE_DIR}/main/tetmagVersion.h.in ${PROJECT_SOURCE_DIR}/main/tetmagVersion.h)

if(CMAKE_BUILD_TYPE MATCHES Release)
  #string(APPEND CMAKE_CXX_FLAGS_RELEASE " -mavx " ) # not compatible with CUDA
  #string(APPEND CMAKE_CXX_FLAGS_RELEASE " -march=native") # not compatible with CUDA
  string(APPEND CMAKE_CXX_FLAGS_RELEASE " -fno-math-errno -Wno-unused -Wall -DNDEBUG")
  message(STATUS "Build flags: ${CMAKE_CXX_FLAGS_RELEASE} ")
endif()

if(CMAKE_BUILD_TYPE MATCHES Debug)
  string(APPEND CMAKE_CXX_FLAGS_DEBUG " -ggdb3")
  message(STATUS "Build flags: ${CMAKE_CXX_FLAGS_DEBUG} ")
endif()

find_package(Eigen3 REQUIRED NO_MODULE)
find_package(Boost CONFIG QUIET COMPONENTS filesystem program_options thread OPTIONAL_COMPONENTS system)
if(NOT Boost_FOUND)
  find_package(Boost REQUIRED COMPONENTS filesystem program_options thread OPTIONAL_COMPONENTS system)
  if(Boost_FOUND AND Boost_VERSION_STRING VERSION_LESS "1.70")
    if(POLICY CMP0167)
      cmake_policy(SET CMP0167 NEW)
    endif()
  endif()
endif()
if(NOT Boost_FOUND)
  message(FATAL_ERROR "Boost not found. Install the Boost library with components: filesystem, program_options, thread.")
endif()
find_package(Threads REQUIRED)
find_package(VTK REQUIRED COMPONENTS CommonCore CommonDataModel IOXML IOLegacy)
find_library(NCDF netcdf REQUIRED)
add_definitions(-DUSE_NETCDF)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${LIBDIR})

set(PREPROC ${PROJECT_SOURCE_DIR}/preproc/)
set(UTILS ${PROJECT_SOURCE_DIR}/utils/)
set(IO ${PROJECT_SOURCE_DIR}/io/)
set(SPECS ${PROJECT_SOURCE_DIR}/specs/)
set(MAIN ${PROJECT_SOURCE_DIR}/main/)
set(GPU ${PROJECT_SOURCE_DIR}/gpu/)

if(USE_CUDA)
  include(${GPU}/CMakeLists.txt)
  add_library(cuda_cccl_includes INTERFACE)
  target_include_directories(cuda_cccl_includes INTERFACE ${CUDAToolkit_INCLUDE_DIRS})

  foreach(root IN LISTS CUDAToolkit_INCLUDE_DIRS)
    if(IS_DIRECTORY "${root}/cccl")
      target_include_directories(cuda_cccl_includes INTERFACE "${root}/cccl")
    endif()
  endforeach()
endif(USE_CUDA)

include(${EXTDIR}/CMakeLists.txt) # install external dependencies

include_directories(${PREPROC} ${UTILS} ${IO} ${SPECS} ${MAIN} ${GPU} ${EIGEN3_INCLUDE_DIR})

set(PREPROC_FILES
  ${PREPROC}BEMprocessing.cpp
  ${PREPROC}BEProperties.cpp
  ${PREPROC}Boundary.cpp
  ${PREPROC}FEMprocessing.cpp
  ${PREPROC}Lindholm.cpp
  ${PREPROC}LindholmCWrapper.cpp
)

set(TETMAG_FILES
  ${PREPROC_FILES}
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
  ${MAIN}LLGIntegrator.cpp
  ${MAIN}SolverFactory.cpp
  ${MAIN}TheLLG.cpp
  ${MAIN}TheSimulation.cpp
  ${UTILS}Timer.cpp
  ${IO}VTKReader.cpp
)

if(USE_CUDA)
  set(CVODE_CUDA sundials_nveccuda)
  add_executable(${PROJECT_NAME} ${TETMAG_FILES} ${CUDA_FILES})
  target_link_libraries(${PROJECT_NAME} PRIVATE cuda_cccl_includes)
  target_link_libraries(${PROJECT_NAME} PRIVATE CUDA::cudart CUDA::cusparse ${CVODE_CUDA})

  target_compile_definitions(${PROJECT_NAME} PRIVATE USE_CUDA)

  target_compile_options(${PROJECT_NAME} PRIVATE
    $<$<COMPILE_LANGUAGE:CUDA>:-Xptxas=-O3>
    $<$<COMPILE_LANGUAGE:CUDA>:-Xcompiler=-O3>
    $<$<COMPILE_LANGUAGE:CUDA>:-Wno-deprecated-gpu-targets>
    $<$<COMPILE_LANGUAGE:CUDA>:-w>
  )
  if(CMAKE_CUDA_HOST_COMPILER)
    target_compile_options(${PROJECT_NAME} PRIVATE  $<$<COMPILE_LANGUAGE:CUDA>:-ccbin=${CMAKE_CUDA_HOST_COMPILER}>)
  endif(CMAKE_CUDA_HOST_COMPILER)
else()
  add_executable(${PROJECT_NAME} ${TETMAG_FILES})
endif(USE_CUDA)

add_dependencies(${PROJECT_NAME} amgcl sundials_cvode h2lib)
target_link_libraries(${PROJECT_NAME} PRIVATE Eigen3::Eigen)
target_link_libraries(${PROJECT_NAME} PRIVATE Threads::Threads)
target_link_libraries(${PROJECT_NAME} PRIVATE Boost::filesystem Boost::program_options Boost::thread)
if(TARGET Boost::system)
  target_link_libraries(${PROJECT_NAME} PRIVATE Boost::system)
endif()
target_link_libraries(${PROJECT_NAME} PRIVATE ${H2LIB} mshio ${CVODELIB})
target_link_libraries(${PROJECT_NAME} PRIVATE VTK::CommonCore VTK::CommonDataModel VTK::IOXML VTK::IOLegacy)
target_link_libraries(${PROJECT_NAME} PRIVATE ${NCDF})

if(USE_OpenMP)
  target_link_libraries(${PROJECT_NAME} PRIVATE OpenMP::OpenMP_CXX)
  set(CVODE_OMP sundials_nvecopenmp)
  add_dependencies(${PROJECT_NAME} ${CVODE_OMP})
  target_link_libraries(${PROJECT_NAME} PRIVATE ${CVODE_OMP})
endif(USE_OpenMP)

set_target_properties(${PROJECT_NAME} PROPERTIES LINKER_LANGUAGE CXX)

target_include_directories(${PROJECT_NAME} PUBLIC ${EIGEN3_INCLUDE_DIR})
target_include_directories(${PROJECT_NAME} PUBLIC $<BUILD_INTERFACE:${EXTDIR}/h2lib/Library>)
target_include_directories(${PROJECT_NAME} PUBLIC $<BUILD_INTERFACE:${amgcl_SOURCE_DIR}>)


if(APPLE)
  set_target_properties(${PROJECT_NAME} PROPERTIES
    INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/${INSTALL_LIB_SUBDIR};/usr/local/lib/"
    BUILD_WITH_INSTALL_RPATH TRUE
  )
else()
  set_target_properties(${PROJECT_NAME} PROPERTIES
    INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/${INSTALL_LIB_SUBDIR}"
    BUILD_WITH_INSTALL_RPATH FALSE
  )
endif()

# Installation
install(TARGETS ${PROJECT_NAME}
  RUNTIME DESTINATION ${INSTALL_BIN_DIR}
)

# Manifest
install(CODE "file(WRITE ${CMAKE_INSTALL_PREFIX}/share/tetmag/install_manifest.txt \"\")")
install(CODE "file(APPEND ${CMAKE_INSTALL_PREFIX}/share/tetmag/install_manifest.txt \"${CMAKE_INSTALL_PREFIX}/${INSTALL_BIN_DIR}/tetmag\\n\")")
install(CODE "file(APPEND ${CMAKE_INSTALL_PREFIX}/share/tetmag/install_manifest.txt \"${CMAKE_INSTALL_PREFIX}/${INSTALL_LIB_SUBDIR}/libh2.a\\n\")")
install(CODE "file(APPEND ${CMAKE_INSTALL_PREFIX}/share/tetmag/install_manifest.txt \"${CMAKE_INSTALL_PREFIX}/${INSTALL_LIB_SUBDIR}/libsundials_cvode.a\\n\")")
if(USE_CUDA)
  install(CODE "file(APPEND ${CMAKE_INSTALL_PREFIX}/share/tetmag/install_manifest.txt \"${CMAKE_INSTALL_PREFIX}/${INSTALL_LIB_SUBDIR}/libsundials_nveccuda.a\\n\")")
endif()
if(USE_OpenMP)
  install(CODE "file(APPEND ${CMAKE_INSTALL_PREFIX}/share/tetmag/install_manifest.txt \"${CMAKE_INSTALL_PREFIX}/${INSTALL_LIB_SUBDIR}/libsundials_nvecopenmp.a\\n\")")
endif()
install(CODE "file(APPEND ${CMAKE_INSTALL_PREFIX}/share/tetmag/install_manifest.txt \"\${HEADER_FILES}\")")

# Uninstall
if(NOT TARGET uninstall)
  add_custom_target(uninstall
    COMMAND ${CMAKE_COMMAND} -E echo "Uninstalling tetmag..."
    COMMAND sh -c "xargs rm -f < ${CMAKE_INSTALL_PREFIX}/share/tetmag/install_manifest.txt"
    COMMAND ${CMAKE_COMMAND} -E remove ${CMAKE_INSTALL_PREFIX}/share/tetmag/install_manifest.txt
    COMMAND sh -c "find ${CMAKE_INSTALL_PREFIX}/${INSTALL_LIB_SUBDIR} -type d -empty -delete"
    COMMAND sh -c "find ${CMAKE_INSTALL_PREFIX}/share/tetmag -type d -empty -delete"
    COMMAND ${CMAKE_COMMAND} -E echo "Uninstallation complete."
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
    COMMENT "Uninstalling tetmag files."
  )
endif()
