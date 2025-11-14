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

#include <algorithm>
#include <cmath>

#include "dragdisk.h"
#include "grid.h"
#include "fields.h"
#include "stats.h"
#include "tools.h"

#include "dragdisk_kernels.cuh"
#include "cuda_launcher.h"
#include "cuda_tiling.h"

namespace
{
    template<typename TF>
    struct DiskRegion
    {
        int k_center = 0;
        int i_center = 0;
        int j_center = 0;
        int istart = 0;
        int iend = 0;
        int jstart = 0;
        int jend = 0;
        TF radius_cells_sq = TF(0);
    };

    template<typename TF>
    DiskRegion<TF> build_disk_region(
            const Grid_data<TF>& gd,
            TF height,
            TF diameter,
            TF xc,
            TF yc)
    {
        DiskRegion<TF> region;

        auto nearest_k = [&](TF target_z) {
            int k_best = gd.kstart;
            TF best = std::abs(gd.z[k_best] - target_z);
            for (int k = gd.kstart + 1; k < gd.kend; ++k)
            {
                TF diff = std::abs(gd.z[k] - target_z);
                if (diff < best)
                {
                    best = diff;
                    k_best = k;
                }
            }
            return k_best;
        };
        region.k_center = nearest_k(height);

        const TF x0 = gd.x[gd.istart];
        const TF y0 = gd.y[gd.jstart];
        const TF inv_dx = TF(1) / gd.dx;
        const TF inv_dy = TF(1) / gd.dy;

        auto clamp_index = [](int val, int low, int high_exclusive)
        {
            if (val < low)
                return low;
            if (val >= high_exclusive)
                return high_exclusive - 1;
            return val;
        };

        const TF offset_i = (xc - x0) * inv_dx;
        const TF offset_j = (yc - y0) * inv_dy;

        region.i_center = clamp_index(
                gd.istart + static_cast<int>(std::lround(offset_i)),
                gd.istart, gd.iend);
        region.j_center = clamp_index(
                gd.jstart + static_cast<int>(std::lround(offset_j)),
                gd.jstart, gd.jend);

        const TF radius_cells = TF(0.5) * diameter / std::sqrt(gd.dx * gd.dy);
        region.radius_cells_sq = radius_cells * radius_cells;

        const int r_int  = static_cast<int>(std::ceil(radius_cells));
        region.istart = std::max(gd.istart, region.i_center - r_int);
        region.iend   = std::min(gd.iend,   region.i_center + r_int + 1);
        region.jstart = std::max(gd.jstart, region.j_center - r_int);
        region.jend   = std::min(gd.jend,   region.j_center + r_int + 1);

        return region;
    }
}

#ifdef USECUDA
template<typename TF>
void DragDisk<TF>::prepare_device()
{
    // No GPU-side allocations are required yet, but keep the hook for future use.
}

template<typename TF>
void DragDisk<TF>::clear_device()
{
}

template<typename TF>
void DragDisk<TF>::exec(Stats<TF>& stats, double)
{
    if (!sw_disk)
        return;

    const auto& gd = grid.get_grid_data();

    const TF diameter_tf = static_cast<TF>(diameter);
    const auto region = build_disk_region(gd, height, diameter_tf, xc, yc);

    k_center = region.k_center;
    radius   = std::sqrt(region.radius_cells_sq);

    Grid_layout grid_layout = {
        region.istart, region.iend,
        region.jstart, region.jend,
        region.k_center, region.k_center + 1,
        gd.istride,
        gd.jstride,
        gd.kstride};

    launch_grid_kernel<DragDisk_kernels::dragdisk_u_g<TF>>(
            grid_layout,
            fields.mt.at("u")->fld_g.view(),
            fields.mp.at("u")->fld_g,
            gd.utrans,
            cd,
            region.i_center,
            region.j_center,
            region.radius_cells_sq);

    cuda_check_error();
    cudaDeviceSynchronize();

    stats.calc_tend(*fields.mt.at("u"), "dragdisk");
}
#endif


#ifdef FLOAT_SINGLE
template class DragDisk<float>;
#else
template class DragDisk<double>;
#endif
