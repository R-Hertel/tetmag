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
################################################################################
# Install H2Lib

if(USE_OpenMP)
  set(H2CONFIG ${CMAKE_COMMAND} -E copy_if_different
    ${EXTDIR}/h2config/options.inc.omp
    ${EXTDIR}/h2lib/options.inc.default
  )
else()
  set(H2CONFIG ${CMAKE_COMMAND} -E copy_if_different
    ${EXTDIR}/h2config/options.inc.serial
    ${EXTDIR}/h2lib/options.inc.default
  )
endif()


if(APPLE)
  if(USE_OpenMP)
    set(H2_CONFIG_FILE ${EXTDIR}/h2config/system.inc.macos)
  else()
    set(H2_CONFIG_FILE ${EXTDIR}/h2lib/system.inc.macos)
  endif()

  set(H2CONFIG ${H2CONFIG} && ${CMAKE_COMMAND} -E copy
    ${H2_CONFIG_FILE}
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
  BUILD_COMMAND ${CMAKE_MAKE_PROGRAM} libh2.a "CC=${CMAKE_C_COMPILER} -std=gnu17"
  BUILD_IN_SOURCE true
  INSTALL_COMMAND cmake -E echo "Skipping install step of h2Lib."
)

add_custom_command(
  TARGET h2lib
  POST_BUILD
  COMMAND ${CMAKE_COMMAND} -E copy ${EXTDIR}/h2lib/libh2.a ${LIBDIR}/
)

add_library(h2 SHARED IMPORTED)
set_target_properties(h2 PROPERTIES IMPORTED_LOCATION ${LIBDIR}/libh2.a)
set(H2LIB h2)
install(FILES ${LIBDIR}/libh2.a DESTINATION ${INSTALL_LIB_SUBDIR})

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

FetchContent_MakeAvailable(amgcl)

###############################################################################
# Install CVODE

include(FetchContent)
FetchContent_Declare(
  cvode
  GIT_REPOSITORY https://github.com/LLNL/sundials.git
  GIT_TAG v6.7.0
  GIT_SHALLOW true
  GIT_PROGRESS true
)


set(SUNDIALS_ENABLE_LAPACK OFF CACHE BOOL "" FORCE)
set(ENABLE_LAPACK OFF CACHE BOOL "" FORCE)
set(ENABLE_SUNLINSOL_DENSE OFF CACHE BOOL "" FORCE)
set(ENABLE_SUNLINSOL_BAND OFF CACHE BOOL "" FORCE)
set(ENABLE_SUNMATRIX_DENSE OFF CACHE BOOL "" FORCE)
set(ENABLE_SUNMATRIX_BAND OFF CACHE BOOL "" FORCE)


set(BUILD_CVODE ON CACHE BOOL "" FORCE)

if(USE_CUDA)
  set(EXAMPLES_ENABLE_CUDA OFF CACHE BOOL "")
  set(CMAKE_CUDA_ARCHITECTURES ${GPU_ARCHITECTURES})
  set(ENABLE_CUDA ON CACHE BOOL "")
  set(BUILD_NVECTOR_CUDA ON CACHE BOOL "")
endif()

if(USE_OpenMP)
  set(ENABLE_OPENMP ON CACHE BOOL "")
  set(BUILD_NVECTOR_OPENMP ON CACHE BOOL "")
endif()

set(BUILD_ARKODE OFF CACHE BOOL "")
set(BUILD_TESTING OFF CACHE BOOL "")
set(BUILD_CVODES OFF CACHE BOOL "")
set(BUILD_KINSOL OFF CACHE BOOL "")
set(BUILD_IDA OFF CACHE BOOL "")
set(BUILD_IDAS OFF CACHE BOOL "")
set(EXAMPLES_ENABLE_C OFF CACHE BOOL "")
set(EXAMPLES_INSTALL OFF CACHE BOOL "")

FetchContent_MakeAvailable(cvode)
add_compile_options(-Wno-reorder)
set(CVODELIB sundials_cvode)

install(FILES ${LIBDIR}/libsundials_cvode.a
  DESTINATION ${INSTALL_LIB_SUBDIR}
)
if(USE_CUDA)
  install(FILES ${LIBDIR}/libsundials_nveccuda.a
    DESTINATION ${INSTALL_LIB_SUBDIR}
  )
endif()
if(USE_OpenMP)
  install(FILES ${LIBDIR}/libsundials_nvecopenmp.a
    DESTINATION ${INSTALL_LIB_SUBDIR}
  )
endif()

###############################################################################
# Install MshIO

FetchContent_Declare(mshio
  GIT_REPOSITORY https://github.com/qnzhou/MshIO.git
  GIT_TAG main
  GIT_SHALLOW true
  GIT_PROGRESS true
)

FetchContent_MakeAvailable(mshio)

################################################################################

