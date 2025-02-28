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
if (TETMAG_FIND_CUDA_ARCH)
  include(FindCUDA/select_compute_arch)
  CUDA_DETECT_INSTALLED_GPUS(INSTALLED_GPU_STRING1)

  string(STRIP "${INSTALLED_GPU_STRING1}" INSTALLED_GPU_STRING2)
  string(REPLACE " " ";" INSTALLED_GPU_STRING3 "${INSTALLED_GPU_STRING2}")
  string(REPLACE "." "" INSTALLED_GPU_ARCHS "${INSTALLED_GPU_STRING3}")
  LIST(APPEND GPU_ARCHITECTURES ${INSTALLED_GPU_ARCHS})
  string(TOLOWER ${INSTALLED_GPU_STRING1} ptx_check)

  string(FIND ${ptx_check} "ptx" ptx_pos) # this string appears in select_compute_arch if no GPU is found
  if (NOT ptx_pos EQUAL -1) # FIND returns a value of -1 if the string is not found.
    message("ERROR: No CUDA capable device found.")
    set (USE_CUDA OFF)
  endif() 
   message(STATUS "Detected GPU architectures: ${INSTALLED_GPU_STRING2}")
else(TETMAG_FIND_CUDA_ARCH)
   set (GPU_ARCHITECTURES ${TETMAG_SET_CUDA_ARCHS})
endif(TETMAG_FIND_CUDA_ARCH)


foreach(ARCH ${GPU_ARCHITECTURES})		
  list(APPEND CUDA_NVCC_FLAGS "-gencode arch=compute_${ARCH},code=sm_${ARCH} ")
endforeach() 

  add_definitions(-DUSE_CUDA)
  set (NVCC_FLAGS "-DUSE_CUDA ; -Xptxas -O3; -Xcompiler -O3; -ccbin ${CMAKE_CUDA_HOST_COMPILER}; ")
  #string(APPEND NVCC_FLAGS "-Xcudafe --diag_suppress=esa_on_defaulted_function_ignored; -Wno-deprecated-gpu-targets; ")
  string(APPEND NVCC_FLAGS " -Wno-deprecated-gpu-targets; ")
  find_package(CUDAToolkit REQUIRED)
  list(APPEND CUDA_NVCC_FLAGS ${NVCC_FLAGS})
 
  set(CMAKE_CUDA_STANDARD 11)

  list(APPEND CUDA_NVCC_FLAGS " -w;")
  message(STATUS "NVCC_FLAGS: " ${CUDA_NVCC_FLAGS}) 

  set (CUDA_FILES
    ${GPU}EffFieldGPU.cu  
    ${GPU}LLGWrapper.cu
    ${GPU}GpuPotential.cu  
    ${GPU}GPUSolver.cu
    ${GPU}SetDevice.cu
    ${GPU}SpMatCUDA.cu
   )
