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

if (TETMAG_FIND_CUDA_ARCH)
  include(FindCUDA/select_compute_arch)
  CUDA_DETECT_INSTALLED_GPUS(INSTALLED_GPU_STRING1)

  string(STRIP "${INSTALLED_GPU_STRING1}" INSTALLED_GPU_STRING2)
  string(REPLACE " " ";" INSTALLED_GPU_STRING3 "${INSTALLED_GPU_STRING2}")
  string(REPLACE "." "" INSTALLED_GPU_ARCHS "${INSTALLED_GPU_STRING3}")
  set(GPU_ARCHITECTURES ${INSTALLED_GPU_ARCHS})

  string(TOLOWER "${INSTALLED_GPU_STRING1}" ptx_check)
  string(FIND "${ptx_check}" "ptx" ptx_pos)
  if (NOT ptx_pos EQUAL -1)
    message(SEND_ERROR "ERROR: No CUDA capable device found.")
    set(USE_CUDA OFF)
  endif()

  message(STATUS "Detected GPU architectures: ${INSTALLED_GPU_STRING2}")
else()
  set(GPU_ARCHITECTURES ${TETMAG_SET_CUDA_ARCHS})
endif()

set(CMAKE_CUDA_ARCHITECTURES "${GPU_ARCHITECTURES}")

set(CUDA_FILES
  ${GPU}EffFieldGPU.cu
  ${GPU}TheLLG_GPU.cu
  ${GPU}DemagGPU.cu
  ${GPU}GPUSolver.cu
  ${GPU}SetDevice.cu
  ${GPU}SpMatCUDA.cu
)
