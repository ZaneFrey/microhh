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
 *
 * MicroHH is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with MicroHH.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "dragdisk.h"

#ifdef USECUDA
#include "tools.h"

namespace {
template<typename TF>
__global__ void dragdisk_kernel(TF* __restrict__ u, const int* __restrict__ idx, int n, TF cd)
{
    int t = blockIdx.x * blockDim.x + threadIdx.x;
    if (t < n)
    {
        int g = idx[t];
        u[g] -= cd * u[g];
    }
}
}

template<typename TF>
void dragdisk_apply_cuda(cuda_vector<TF>& u, const cuda_vector<int>& indices, TF cd)
{
    int n = indices.size();
    if (n == 0) return;
    int block = 128;
    int grid = (n + block - 1) / block;
    dragdisk_kernel<<<grid, block>>>(u.data(), indices.data(), n, cd);
    cuda_check_error();
}

template void dragdisk_apply_cuda<float>(cuda_vector<float>&, const cuda_vector<int>&, float);
template void dragdisk_apply_cuda<double>(cuda_vector<double>&, const cuda_vector<int>&, double);
#endif