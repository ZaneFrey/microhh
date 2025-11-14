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
#include <vector>

#include "dragdisk.h"
#include "grid.h"
#include "fields.h"
#include "field3d_operators.h"
#include "fast_math.h"
#include "input.h"
#include "master.h"
#include "netcdf_interface.h"
#include "stats.h"

namespace
{
    namespace fm = Fast_math;

    template<typename TF>
    void find_disk_points(
        const Grid<TF>& grid,
        TF height,
        TF diameter,
        TF xc,
        TF yc,
        int& k_center,
        int& i_center,
        int& j_center,
        TF& radius_cells,
        std::vector<int>& indices_out)
    {
        const auto& gd = grid.get_grid_data();
        const int jj = gd.jstride;
        const int kk = gd.kstride;

        // Snap the requested height to the nearest full level
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
        k_center = nearest_k(height);

        // Locate the closest i/j centers
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

        i_center = clamp_index(
                gd.istart + static_cast<int>(std::lround(offset_i)),
                gd.istart, gd.iend);
        j_center = clamp_index(
                gd.jstart + static_cast<int>(std::lround(offset_j)),
                gd.jstart, gd.jend);

        // Convert physical diameter to a radius in index space
        const TF spacing = TF(0.5) * diameter / std::sqrt(gd.dx * gd.dy);
        radius_cells = TF(0.5) * diameter / spacing;

        // Scan the bounding box around the disk and keep points inside the circle
        const int r_int  = static_cast<int>(std::ceil(radius_cells));
        const int istart = std::max(gd.istart, i_center - r_int);
        const int iend   = std::min(gd.iend,   i_center + r_int + 1);
        const int jstart = std::max(gd.jstart, j_center - r_int);
        const int jend   = std::min(gd.jend,   j_center + r_int + 1);

        indices_out.clear();
        for (int j = jstart; j < jend; ++j)
            for (int i = istart; i < iend; ++i)
            {
                const TF di = TF(i - i_center);
                const TF dj = TF(j - j_center);
                if (di * di + dj * dj <= radius_cells * radius_cells)
                {
                    int idx = i + j * jj + k_center * kk;
                    indices_out.push_back(idx);
                }
            }

    }

    template<typename TF>
    void dragdisk_u(
            TF* const restrict ut,
            const TF* const restrict u,
            const TF* const restrict v,
            const TF* const restrict w,
            const TF utrans,
            const TF vtrans,
            const TF cd,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int kstart, const int kend,
            const int jstride, const int kstride)    
    {
        const int ii = 1;
        const int jj = jstride;
        const int kk = kstride;
        for (int k=kstart; k<kend; ++k)
            for (int j=jstart; j<jend; ++j)
                for (int i=istart; i<iend; ++i)
                {
                    const int ijk = i + j*jj + k*kk;
                    const TF u_on_u = u[ijk] + utrans;
                    const TF ftau = -cd * std::abs(u_on_u);
                    
                    ut[ijk] += ftau * u_on_u;                
                }
    }
} 

template<typename TF>
DragDisk<TF>::DragDisk(
        Master& masterin, Grid<TF>& gridin, Fields<TF>& fieldsin, Input& input) :
        master(masterin), grid(gridin), fields(fieldsin), 
        field3d_operators(masterin, gridin, fieldsin), field3d_io(masterin, gridin)
{
    sw_disk = input.get_item<bool>("dragdisk", "swdragdisk", "", false);

    if (sw_disk)
    {
        diameter = input.get_item<int>("dragdisk", "diameter", "", TF(100));
        height   = input.get_item<TF>("dragdisk", "height", "", TF(100));
        cd       = input.get_item<TF>("dragdisk", "cd", "");
        xc       = input.get_item<TF>("dragdisk", "xc", "", TF(100));
        yc       = input.get_item<TF>("dragdisk", "yc", "", TF(100));
    }
}

template<typename TF>
DragDisk<TF>::~DragDisk()
{
}

template <typename TF>
void DragDisk<TF>::init(
        Input& inputin)
{
    if (!sw_disk)
        return;

    auto& gd = grid.get_grid_data();

}

template<typename TF>
void DragDisk<TF>::create(
        Input& inputin, Stats<TF>& stats)
{
    if (!sw_disk)
        return;

    auto& gd = grid.get_grid_data();

    disk_indices.clear();
    find_disk_points(grid, height, diameter, xc, yc, k_center, i_center, j_center, radius, disk_indices);
}

#ifndef USECUDA
template<typename TF>
void DragDisk<TF>::exec()
{
    if (!sw_disk)
        return;

    auto& gd = grid.get_grid_data();

    dragdisk_u<TF>(
        fields.mt.at("u")->fld.data(),
        fields.mp.at("u")->fld.data(),
        fields.mp.at("v")->fld.data(),
        fields.mp.at("w")->fld.data(),
        gd.utrans,
        gd.vtrans,
        cd,
        gd.istart, gd.iend,
        gd.jstart, gd.jend,
        k_center, k_center + 1,
        gd.jstride, gd.kstride);
}
#endif


#ifdef FLOAT_SINGLE
template class DragDisk<float>;
#else
template class DragDisk<double>;
#endif
