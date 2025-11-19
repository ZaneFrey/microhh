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
    
}

#ifdef USECUDA
template<typename TF>
void DragDisk<TF>::prepare_device()
{
}

template<typename TF>
void DragDisk<TF>::clear_device()
{
}

template<typename TF>
void DragDisk<TF>::exec(Stats<TF>& stats)
{
    if (!sw_disk)
        return;

    // If the disk does not intersect this MPI rank in y-z, CPU-side
    // create() left disk_indices empty; then there is nothing to do here.
    if (disk_indices.empty())
        return;

    const auto& gd = grid.get_grid_data();

    // Reconstruct the vertical index closest to the disk center height.
    int k_center = gd.kstart;
    TF best = std::abs(gd.z[k_center] - height);
    for (int k = gd.kstart + 1; k < gd.kend; ++k)
    {
        const TF diff = std::abs(gd.z[k] - height);
        if (diff < best)
        {
            best = diff;
            k_center = k;
        }
    }

    // Compute effective radius following the CPU implementation:
    // use a reduced radius so that only grid cells fully inside the disk
    // are included.
    const TF radius = TF(0.5) * diameter;

    TF dz_center;
    if (k_center > gd.kstart && k_center < gd.kend - 1)
        dz_center = TF(0.5) * (gd.z[k_center + 1] - gd.z[k_center - 1]);
    else if (k_center < gd.kend - 1)
        dz_center = gd.z[k_center + 1] - gd.z[k_center];
    else
        dz_center = gd.z[k_center] - gd.z[k_center - 1];

    const TF cell_diag   = std::sqrt(gd.dy * gd.dy + dz_center * dz_center);
    const TF eff_radius  = std::max(TF(0), radius - TF(0.5) * cell_diag);
    const TF eff_radius2 = eff_radius * eff_radius;

    // Restrict the launch to the rotor i-plane at i_center.
    Grid_layout grid_layout_disk = {
            i_center, i_center + 1,
            gd.jstart, gd.jend,
            gd.kstart, gd.kend,
            gd.istride,
            gd.jstride,
            gd.kstride};

    launch_grid_kernel<DragDisk_kernels::dragdisk_u_g<TF>>(
            grid_layout_disk,
            fields.mt.at("u")->fld_g.view(),
            fields.mp.at("u")->fld_g,
            gd.utrans,
            cd,
            gd.y_g,
            gd.z_g,
            yc,
            height,
            eff_radius2);

    cudaDeviceSynchronize();
    stats.calc_tend(*fields.mt.at("u"), tend_name);
}
#endif


#ifdef FLOAT_SINGLE
template class DragDisk<float>;
#else
template class DragDisk<double>;
#endif
