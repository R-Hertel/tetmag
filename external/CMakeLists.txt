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
################################################################################
# Install H2Lib


if (USE_OpenMP) 
  set (H2CONFIG ${CMAKE_COMMAND} -E copy_if_different
    ${EXTDIR}/h2config/options.inc.omp
    ${EXTDIR}/h2lib/options.inc.default
    )
Else()
  set (H2CONFIG ${CMAKE_COMMAND} -E copy_if_different
    ${EXTDIR}/h2config/options.inc.serial
    ${EXTDIR}/h2lib/options.inc.default
    )
endif()

if (APPLE)
  set(H2CONFIG ${H2CONFIG} && ${CMAKE_COMMAND} -E copy 
    ${EXTDIR}/h2lib/system.inc.macos 
    ${EXTDIR}/h2lib/system.inc.linux
  )
endif()

include(ExternalProject)
ExternalProject_Add(h2lib
   GIT_REPOSITORY https://github.com/H2Lib/H2Lib.git
   GIT_TAG collocation
   UPDATE_DISCONNECTED true
   SOURCE_DIR ${EXTDIR}/h2lib
   CONFIGURE_COMMAND ${H2CONFIG}     
   BUILD_IN_SOURCE true
   INSTALL_COMMAND cmake -E echo "Skipping install step of h2Lib."
)

add_custom_command(
    TARGET h2lib
    POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy ${EXTDIR}/h2lib/libh2.a  ${LIBDIR}/
)

add_library(h2 SHARED IMPORTED)
set_target_properties(h2 PROPERTIES IMPORTED_LOCATION ${LIBDIR}/libh2.a)
set(H2LIB h2)

###############################################################################
# Istall AMGCL

include(FetchContent)
FetchContent_Declare(amgcl
  GIT_REPOSITORY https://github.com/ddemidov/amgcl.git
  GIT_TAG master
  GIT_SHALLOW true
  GIT_PROGRESS true    
  )

set(AMGCL_BUILD_EXAMPLES OFF CACHE BOOL "")
set(AMGCL_BUILD_TESTS OFF CACHE BOOL "")
FetchContent_GetProperties(amgcl)

if(NOT amgcl_POPULATED)
  FetchContent_Populate(amgcl)
  add_subdirectory(${amgcl_SOURCE_DIR} ${EXTDIR}/amgcl)
endif()
FetchContent_MakeAvailable(amgcl)

###############################################################################
# Install CVODE

include(FetchContent)
if(TETMAG_USE_CVODE_5)   
  FetchContent_Declare(cvode
    GIT_REPOSITORY https://github.com/LLNL/sundials.git
    GIT_TAG v5.1.0
    GIT_SHALLOW true
    GIT_PROGRESS true
    )
else()
    FetchContent_Declare(cvode
    GIT_REPOSITORY https://github.com/LLNL/sundials.git
    GIT_TAG v6.7.0
#    GIT_TAG main
    GIT_SHALLOW true
    GIT_PROGRESS true
    )
endif()

set(BUILD_CVODE ON)

if (USE_CUDA)
  set(EXAMPLES_ENABLE_CUDA OFF CACHE BOOL "")
  set(CMAKE_CUDA_ARCHITECTURES ${GPU_ARCHITECTURES})  
  set(ENABLE_CUDA ON CACHE BOOL "")
  set(BUILD_NVECTOR_CUDA ON CACHE BOOL "")
  set(CMAKE_CUDA_STANDARD 14)
endif()
  
if (USE_OpenMP)
  set(ENABLE_OPENMP ON CACHE BOOL "")
  set(BUILD_NVECTOR_OPENMP ON CACHE BOOL "")
endif()

set(BUILD_ARKODE OFF CACHE BOOL "")
set(BUILD_TESTING OFF CACHE BOOL "")
set(BUILD_CVODES OFF CACHE BOOL "")
set(BUILD_KINSOL OFF CACHE BOOL "")
set(BUILD_IDA OFF CACHE BOOL "")
set(BUILD_IDAS OFF CACHE BOOL "" )
set(EXAMPLES_ENABLE_C OFF CACHE BOOL "")
set(EXAMPLES_INSTALL OFF CACHE BOOL "")

FetchContent_GetProperties(cvode)

if(NOT cvode_POPULATED)
  FetchContent_Populate(cvode)
  add_subdirectory(${cvode_SOURCE_DIR} ${EXTDIR}/sundials)
endif()
FetchContent_MakeAvailable(cvode)
set(CVODELIB sundials_cvode)

###############################################################################
# Intall GMSH


if(WIN32)
  if (CMAKE_SIZEOF_VOID_P EQUAL 8)
    set(GMSH_DOWNLOAD_LINK "https://gmsh.info/bin/Windows/gmsh-git-Windows64-sdk.zip")
  else()
    message(FATAL_ERROR "### ERROR: The 32-bit version of Windows is not supported.")
  endif()
endif(WIN32)

if(APPLE)
  if(CMAKE_SYSTEM_PROCESSOR STREQUAL "arm64")
   set(GMSH_DOWNLOAD_LINK "https://gmsh.info/bin/macOS/gmsh-git-MacOSARM-sdk.tgz")
  else()
   set(GMSH_DOWNLOAD_LINK "https://gmsh.info/bin/macOS/gmsh-git-MacOSX-sdk.tgz")
  endif()
endif(APPLE)

if(UNIX AND NOT APPLE)
  set(GMSH_DOWNLOAD_LINK "https://gmsh.info/bin/Linux/gmsh-git-Linux64-sdk.tgz")
endif(UNIX AND NOT APPLE)

include(ExternalProject)
ExternalProject_Add(gmsh
  URL ${GMSH_DOWNLOAD_LINK}
  UPDATE_DISCONNECTED true
  SOURCE_DIR ${EXTDIR}/gmsh
  CONFIGURE_COMMAND ""
  BUILD_COMMAND ""
  INSTALL_COMMAND ""
  )

add_custom_command(
    TARGET gmsh
    POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy ${EXTDIR}/gmsh/lib/libgmsh*  ${LIBDIR}/
)

add_library(gmsh_lib SHARED IMPORTED)
if (APPLE)
   set_target_properties(gmsh_lib PROPERTIES IMPORTED_LOCATION ${LIBDIR}/libgmsh.dylib)
endif(APPLE)
if (UNIX AND NOT APPLE)
   set_target_properties(gmsh_lib PROPERTIES IMPORTED_LOCATION ${LIBDIR}/libgmsh.so)
endif(UNIX AND NOT APPLE)
set(GMSHLIB gmsh_lib)
################################################################################

