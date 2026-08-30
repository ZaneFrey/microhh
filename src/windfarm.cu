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

#include <cuda_runtime.h>

#include "fields.h"
#include "grid.h"
#include "master.h"
#include "tools.h"
#include "windfarm.h"
#include "windfarm_kernels.cuh"

namespace
{
    template<typename TF>
    __device__ TF interpolate_g(
            const TF* const __restrict__ field,
            const Windfarm_interpolation<TF>& interpolation)
    {
        TF value = TF(0);
        #pragma unroll
        for (int n=0; n<8; ++n)
            value += interpolation.weight[n]*field[interpolation.index[n]];
        return value;
    }

    template<typename TF>
    __global__ void sample_g(
            TF* const __restrict__ sums,
            const Windfarm_sample<TF>* const __restrict__ samples,
            const int sample_count,
            const TF* const __restrict__ u,
            const TF* const __restrict__ v,
            const TF* const __restrict__ w,
            const TF utrans,
            const TF vtrans)
    {
        const int p = blockIdx.x*blockDim.x+threadIdx.x;
        if (p >= sample_count || !samples[p].owner)
            return;
        const Windfarm_sample<TF>& sample = samples[p];
        const TF velocity = (interpolate_g(u, sample.u)+utrans)*sample.nx+
                            (interpolate_g(v, sample.v)+vtrans)*sample.ny+
                            TF(0)*interpolate_g(w, sample.w);
        atomicAdd(&sums[3*sample.turbine], velocity*sample.area);
        atomicAdd(&sums[3*sample.turbine+1], sample.density*sample.area);
        atomicAdd(&sums[3*sample.turbine+2], sample.area);
    }

    template<typename TF>
    __global__ void scatter_g(
            TF* const __restrict__ source,
            const Windfarm_scatter<TF>* const __restrict__ scatter,
            const int scatter_count,
            const Windfarm_force<TF>* const __restrict__ force,
            const int component)
    {
        const int n = blockIdx.x*blockDim.x+threadIdx.x;
        if (n >= scatter_count)
            return;
        const Windfarm_scatter<TF>& entry = scatter[n];
        const TF value = component == 0 ? force[entry.element].x :
                         component == 1 ? force[entry.element].y : force[entry.element].z;
        atomicAdd(&source[entry.index], value*entry.coefficient);
    }

    template<typename TF>
    __global__ void add_source_g(
            TF* const __restrict__ tendency,
            const TF* const __restrict__ source,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int kstart, const int kend,
            const int icells, const int ijcells)
    {
        const int i = blockIdx.x*blockDim.x+threadIdx.x+istart;
        const int j = blockIdx.y*blockDim.y+threadIdx.y+jstart;
        const int k = blockIdx.z+kstart;
        if (i < iend && j < jend && k < kend)
        {
            const int index = i+j*icells+k*ijcells;
            tendency[index] += source[index];
        }
    }
}

template<typename TF>
void launch_windfarm_sample_g(
        TF* const sums, const Windfarm_sample<TF>* const samples, const int sample_count,
        const TF* const u, const TF* const v, const TF* const w, const int turbine_count,
        const TF utrans, const TF vtrans)
{
    cuda_safe_call(cudaMemset(sums, 0, 3*turbine_count*sizeof(TF)));
    if (sample_count > 0)
        sample_g<TF><<<(sample_count+255)/256, 256>>>(sums, samples, sample_count, u, v, w, utrans, vtrans);
    cuda_check_error();
}

template<typename TF>
void launch_windfarm_scatter_g(
        TF* const source, const Windfarm_scatter<TF>* const scatter, const int scatter_count,
        const Windfarm_force<TF>* const force, const int ncells, const int component)
{
    cuda_safe_call(cudaMemset(source, 0, ncells*sizeof(TF)));
    if (scatter_count > 0)
        scatter_g<TF><<<(scatter_count+255)/256, 256>>>(source, scatter, scatter_count, force, component);
    cuda_check_error();
}

template<typename TF>
void launch_windfarm_add_source_g(
        TF* const tendency, const TF* const source,
        const int istart, const int iend, const int jstart, const int jend,
        const int kstart, const int kend, const int icells, const int ijcells)
{
    const dim3 block(16, 16, 1);
    const dim3 grid((iend-istart+15)/16, (jend-jstart+15)/16, kend-kstart);
    add_source_g<TF><<<grid, block>>>(
            tendency, source, istart, iend, jstart, jend, kstart, kend, icells, ijcells);
    cuda_check_error();
}

template<typename TF>
void WindFarm<TF>::prepare_device()
{
    if (!created)
        return;
    samples_g = cuda_vector<Windfarm_sample<TF>>(samples);
    element_forces_g = cuda_vector<Windfarm_force<TF>>(element_forces);
    sums_g.allocate(sums.size());
    for (int component=0; component<3; ++component)
    {
        scatter_g[component] = cuda_vector<Windfarm_scatter<TF>>(scatter[component]);
        source_g[component].allocate(source[component].size());
        cuda_safe_call(cudaMemset(source_g[component].data(), 0, source_g[component].size()*sizeof(TF)));
    }
    device_prepared = true;
}

template<typename TF>
void WindFarm<TF>::clear_device()
{
    samples_g.free();
    element_forces_g.free();
    sums_g.free();
    for (int component=0; component<3; ++component)
    {
        scatter_g[component].free();
        source_g[component].free();
    }
    device_prepared = false;
}

template<typename TF>
void WindFarm<TF>::update_g()
{
    if (!device_prepared)
        throw std::runtime_error("Wind-farm device data have not been prepared");
    const auto& gd = grid.get_grid_data();
    launch_windfarm_sample_g(
            sums_g.data(), samples_g.data(), int(samples_g.size()),
            fields.mp.at("u")->fld_g.data(), fields.mp.at("v")->fld_g.data(),
            fields.mp.at("w")->fld_g.data(), int(turbines.size()), gd.utrans, gd.vtrans);
    cuda_copy(sums_g, sums);
    master.sum(sums.data(), int(sums.size()));
}

template<typename TF>
void WindFarm<TF>::exec_g()
{
    if (!created)
        return;
    const auto& gd = grid.get_grid_data();
    launch_windfarm_add_source_g(
            fields.mt.at("u")->fld_g.data(), source_g[0].data(), gd.istart, gd.iend,
            gd.jstart, gd.jend, gd.kstart, gd.kend, gd.icells, gd.ijcells);
    launch_windfarm_add_source_g(
            fields.mt.at("v")->fld_g.data(), source_g[1].data(), gd.istart, gd.iend,
            gd.jstart, gd.jend, gd.kstart, gd.kend, gd.icells, gd.ijcells);
    launch_windfarm_add_source_g(
            fields.mt.at("w")->fld_g.data(), source_g[2].data(), gd.istart, gd.iend,
            gd.jstart, gd.jend, gd.kstart+1, gd.kend, gd.icells, gd.ijcells);
}

template<typename TF>
void WindFarm<TF>::exec()
{
    exec_g();
}

#ifdef FLOAT_SINGLE
template void launch_windfarm_sample_g<float>(float*, const Windfarm_sample<float>*, int,
        const float*, const float*, const float*, int, float, float);
template void launch_windfarm_scatter_g<float>(float*, const Windfarm_scatter<float>*, int,
        const Windfarm_force<float>*, int, int);
template void launch_windfarm_add_source_g<float>(float*, const float*, int, int, int, int, int, int, int, int);
template class WindFarm<float>;
#else
template void launch_windfarm_sample_g<double>(double*, const Windfarm_sample<double>*, int,
        const double*, const double*, const double*, int, double, double);
template void launch_windfarm_scatter_g<double>(double*, const Windfarm_scatter<double>*, int,
        const Windfarm_force<double>*, int, int);
template void launch_windfarm_add_source_g<double>(double*, const double*, int, int, int, int, int, int, int, int);
template class WindFarm<double>;
#endif
