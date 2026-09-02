#
#  MicroHH
#  Copyright (c) 2011-2024 Chiel van Heerwaarden
#  Copyright (c) 2011-2024 Thijs Heus
#  Copyright (c) 2014-2024 Bart van Stratum
#
#  This file is part of MicroHH
#
#  MicroHH is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  MicroHH is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with MicroHH.  If not, see <http://www.gnu.org/licenses/>.
#
FILE(GLOB sourcefiles "../src/*.cxx" "../src/*.cu")
include_directories("../include" "../rte-rrtmgp-cpp/include" SYSTEM ${INCLUDE_DIRS})

add_library(microhhc STATIC ${sourcefiles})
target_include_directories(microhhc PRIVATE "../include" "../rte-rrtmgp-cpp/include" "../rte-rrtmgp-cpp/include_rt" "../rte-rrtmgp-cpp/include_kernels_cuda" "../rte-rrtmgp-cpp/include_rt_kernels" ${INCLUDE_DIRS})

if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    set_source_files_properties(diff_amd2.cxx PROPERTIES
        COMPILE_OPTIONS "-fno-fast-math;-fno-finite-math-only;-fno-associative-math;-ffp-contract=off")
elseif(CMAKE_CXX_COMPILER_ID MATCHES "Intel|IntelLLVM")
    set_source_files_properties(diff_amd2.cxx PROPERTIES COMPILE_OPTIONS "-fp-model=strict")
else()
    message(FATAL_ERROR "AMD2 has no verified strict-FP policy for C++ compiler ${CMAKE_CXX_COMPILER_ID}")
endif()

if(USECUDA)
    # AMD coefficient/status boundaries require IEEE classification and the
    # exact unfused operation order documented by the model formulation.
    set_source_files_properties(diff_amd2.cu PROPERTIES
        COMPILE_OPTIONS "--fmad=false;--prec-div=true;--prec-sqrt=true;--ftz=false")
    if(USEKERNELLAUNCHER)
        target_link_libraries(microhhc PRIVATE CUDA::curand CUDA::cufft ${LIBS} kernel_launcher)
    else()
        target_link_libraries(microhhc PRIVATE CUDA::curand CUDA::cufft ${LIBS})
    endif()
else()
    target_link_libraries(microhhc PRIVATE ${LIBS})
endif()
