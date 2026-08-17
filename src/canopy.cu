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

#include <stdexcept>

#include "canopy.h"
#include "fields.h"
#include "grid.h"
#include "tools.h"

#include "canopy_kernels.cuh"
#include "cuda_launcher.h"
#include "cuda_tiling.h"

#ifdef USECUDA
template<typename TF>
void Canopy<TF>::prepare_device()
{
    if (!sw_canopy)
        return;

    if (sw_3d_pad)
        throw std::runtime_error("3D plant area density is not yet implemented on GPU.");

    if (pad.empty() || padh.empty())
        throw std::runtime_error("Canopy: prepare_device called before PAD profiles were allocated/loaded.");

    if (pad_g != nullptr)
        cuda_safe_call(cudaFree(pad_g));
    if (padh_g != nullptr)
        cuda_safe_call(cudaFree(padh_g));

    cuda_safe_call(cudaMalloc(reinterpret_cast<void**>(&pad_g), pad.size() * sizeof(TF)));
    cuda_safe_call(cudaMemcpy(pad_g, pad.data(), pad.size() * sizeof(TF), cudaMemcpyHostToDevice));

    cuda_safe_call(cudaMalloc(reinterpret_cast<void**>(&padh_g), padh.size() * sizeof(TF)));
    cuda_safe_call(cudaMemcpy(padh_g, padh.data(), padh.size() * sizeof(TF), cudaMemcpyHostToDevice));
}

template<typename TF>
void Canopy<TF>::clear_device()
{
    if (pad_g != nullptr)
    {
        cuda_safe_call(cudaFree(pad_g));
        pad_g = nullptr;
    }

    if (padh_g != nullptr)
    {
        cuda_safe_call(cudaFree(padh_g));
        padh_g = nullptr;
    }
}

template<typename TF>
void Canopy<TF>::exec()
{
    if (!sw_canopy)
        return;

    if (sw_3d_pad)
        throw std::runtime_error("3D plant area density is not yet implemented on GPU.");

    if (pad_g == nullptr || padh_g == nullptr)
        throw std::runtime_error("Canopy: device PAD not available; did you call prepare_device()?");

    const auto& gd = grid.get_grid_data();

    Grid_layout grid_layout_u = {
            gd.istart, gd.iend,
            gd.jstart, gd.jend,
            gd.kstart, kend_canopy,
            gd.istride,
            gd.jstride,
            gd.kstride};

    launch_grid_kernel<Canopy_kernels::canopy_drag_u_g<TF>>(
            grid_layout_u,
            fields.mt.at("u")->fld_g.view(),
            fields.mp.at("u")->fld_g,
            fields.mp.at("v")->fld_g,
            fields.mp.at("w")->fld_g,
            pad_g,
            gd.utrans,
            gd.vtrans,
            cd);

    launch_grid_kernel<Canopy_kernels::canopy_drag_v_g<TF>>(
            grid_layout_u,
            fields.mt.at("v")->fld_g.view(),
            fields.mp.at("u")->fld_g,
            fields.mp.at("v")->fld_g,
            fields.mp.at("w")->fld_g,
            pad_g,
            gd.utrans,
            gd.vtrans,
            cd);

    const int kstart_w = gd.kstart + 1;
    if (kstart_w < kend_canopy)
    {
        Grid_layout grid_layout_w = {
                gd.istart, gd.iend,
                gd.jstart, gd.jend,
                kstart_w,  kend_canopy,
                gd.istride,
                gd.jstride,
                gd.kstride};

        launch_grid_kernel<Canopy_kernels::canopy_drag_w_g<TF>>(
                grid_layout_w,
                fields.mt.at("w")->fld_g.view(),
                fields.mp.at("u")->fld_g,
                fields.mp.at("v")->fld_g,
                fields.mp.at("w")->fld_g,
                padh_g,
                gd.utrans,
                gd.vtrans,
                cd);
    }
}
#endif

#ifdef FLOAT_SINGLE
template class Canopy<float>;
#else
template class Canopy<double>;
#endif
