/*
 * MicroHH
 * Copyright (c) 2011-2024 Chiel van Heerwaarden
 * Copyright (c) 2011-2024 Thijs Heus
 * Copyright (c) 2014-2024 Bart van Stratum
 *
 * This file is part of MicroHH
 *
 * MicroHH is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.

 * MicroHH is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.

 * You should have received a copy of the GNU General Public License
 * along with MicroHH.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "boundary_cyclic_cuda.cuh"
#include "tools.h"

namespace
{
    template<typename T>
    __global__ void cyclic_x_g(
            T* const data, const int icells, const int jcells, const int kcells,
            const int istart, const int iend)
    {
        const int index = blockIdx.x*blockDim.x+threadIdx.x;
        const int count = istart*jcells*kcells;
        if (index < count)
        {
            const int i = index%istart;
            const int remainder = index/istart;
            const int j = remainder%jcells;
            const int k = remainder/jcells;
            const int stride = icells*jcells;
            data[i+j*icells+k*stride] = data[iend-istart+i+j*icells+k*stride];
            data[iend+i+j*icells+k*stride] = data[istart+i+j*icells+k*stride];
        }
    }

    template<typename T>
    __global__ void cyclic_y_g(
            T* const data, const int icells, const int jcells, const int kcells,
            const int jstart, const int jend)
    {
        const int index = blockIdx.x*blockDim.x+threadIdx.x;
        const int count = jstart*icells*kcells;
        if (index < count)
        {
            const int i = index%icells;
            const int remainder = index/icells;
            const int j = remainder%jstart;
            const int k = remainder/jstart;
            const int stride = icells*jcells;
            if (jend-jstart > 1)
            {
                data[i+j*icells+k*stride] = data[i+(jend-jstart+j)*icells+k*stride];
                data[i+(jend+j)*icells+k*stride] = data[i+(jstart+j)*icells+k*stride];
            }
            else
            {
                data[i+j*icells+k*stride] = data[i+jstart*icells+k*stride];
                data[i+(jend+j)*icells+k*stride] = data[i+jstart*icells+k*stride];
            }
        }
    }
}

template<typename T>
void launch_boundary_cyclic_g(
        T* const data, const int icells, const int jcells, const int kcells,
        const int istart, const int iend, const int jstart, const int jend)
{
    const int x_count = istart*jcells*kcells;
    cyclic_x_g<<<(x_count+255)/256, 256>>>(data, icells, jcells, kcells, istart, iend);
    const int y_count = jstart*icells*kcells;
    cyclic_y_g<<<(y_count+255)/256, 256>>>(data, icells, jcells, kcells, jstart, jend);
    cuda_check_error();
}

template void launch_boundary_cyclic_g<float>(float*, int, int, int, int, int, int, int);
template void launch_boundary_cyclic_g<double>(double*, int, int, int, int, int, int, int);
template void launch_boundary_cyclic_g<unsigned int>(unsigned int*, int, int, int, int, int, int, int);
