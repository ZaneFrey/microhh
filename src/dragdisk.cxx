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
#include <stdexcept>

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
    // TODO: Need to add turbine index into the I/O of the find_disk_points for multiple turbine cases
    // TODO: Need to find fractional area for points not entirely within disk extent
    template<typename TF>
    void find_disk_points(
        const Grid<TF>& grid,
        TF height,
        TF diameter,
        TF xc,
        TF yc,
        int& i_center,
        std::vector<int>& disk_indices)
    {
        const auto& gd = grid.get_grid_data();

        // Ensure the disk lies inside the physical domain.
        const TF radius = TF(0.5) * diameter;
        if (yc - radius < TF(0) || yc + radius > gd.ysize ||
            height - radius < TF(0) || height + radius > gd.zsize)
        {
            throw std::runtime_error("DragDisk: disk boundaries lie outside the global domain");
        }

        // Find the disk center height index.
        int    k_center = gd.kstart;
        TF     k_best   = std::abs(gd.z[k_center] - height);
        for (int k = gd.kstart + 1; k < gd.kend; ++k)
        {
            const TF diff = std::abs(gd.z[k] - height);
            if (diff < k_best)
            {
                k_best   = diff;
                k_center = k;
            }
        }

        // Convert the disk center coordinate (x [m], y[m]) to local grid indices.
        const TF x0 = gd.x[gd.istart];
        const TF y0 = gd.y[gd.jstart];

        const TF fi = (xc - x0) * gd.dxi;
        const TF fj = (yc - y0) * gd.dyi;

        int i_cand = gd.istart + static_cast<int>(std::floor(fi + TF(0.5)));
        int j_center = gd.jstart + static_cast<int>(std::floor(fj + TF(0.5)));

        i_cand   = std::max(gd.istart, std::min(gd.iend - 1, i_cand));
        j_center = std::max(gd.jstart, std::min(gd.jend - 1, j_center));

        // Check whether this rank actually owns the x-location of the disk.
        const TF x_center = gd.x[i_cand];
        if (std::abs(x_center - xc) > TF(0.5) * gd.dx)
        {
            disk_indices.clear();
            return;
        }
        i_center = i_cand;

        // Compute disk extent in indices (used only as a bounding box).
        const TF dy = gd.dy;

        int j_extent = static_cast<int>(std::ceil(radius * gd.dyi));

        TF dz_center;
        if (k_center > gd.kstart && k_center < gd.kend - 1)
            dz_center = TF(0.5) * (gd.z[k_center + 1] - gd.z[k_center - 1]); // interior cell
        else if (k_center < gd.kend - 1)
            dz_center = gd.z[k_center + 1] - gd.z[k_center]; // bottom edge case
        else
            dz_center = gd.z[k_center] - gd.z[k_center - 1]; // top edge case

        int k_extent = static_cast<int>(std::ceil(radius / dz_center));

        const int j_min = j_center - j_extent;
        const int j_max = j_center + j_extent;
        const int k_min = k_center - k_extent;
        const int k_max = k_center + k_extent;

        // Intersect the disk bounding box with the local sub-domain.
        const int j_lo = std::max(j_min, gd.jstart);        // CHECK THIS
        const int j_hi = std::min(j_max, gd.jend - 1);      // CHECK THIS for MPI issue
        const int k_lo = std::max(k_min, gd.kstart);
        const int k_hi = std::min(k_max, gd.kend - 1);

        disk_indices.clear();

        if (j_lo > j_hi || k_lo > k_hi)
            return; // Disk does not intersect this rank in y-z.

        // Extract disk indices on y-z plane that fall entirely within the diameter.
        // Indices are flattened as j * jstride + k * kstride so they can be reused in dragdisk_u.
        const TF cell_diag = std::sqrt(dy * dy + dz_center * dz_center);
        const TF eff_radius = std::max(TF(0), radius - TF(0.5) * cell_diag);
        const TF eff_radius2 = eff_radius * eff_radius;

        const int ny = j_hi - j_lo + 1;
        const int nz = k_hi - k_lo + 1;
        disk_indices.reserve(ny * nz);

        const int jstride = gd.jstride;
        const int kstride = gd.kstride;

        for (int k = k_lo; k <= k_hi; ++k)
        {
            const TF dz = TF(k - k_center) * dz_center;
            const TF dz2 = dz * dz;

            for (int j = j_lo; j <= j_hi; ++j)
            {
                const TF dy_off = TF(j - j_center) * dy;
                const TF dist2 = dy_off * dy_off + dz2;

                if (dist2 <= eff_radius2)
                {
                    const int index = j * jstride + k * kstride;
                    disk_indices.push_back(index);
                }
            }
        }
    }

    template<typename TF>
    void dragdisk_u(
            TF* const restrict ut,
            const TF* const restrict u,
            const TF utrans,
            const TF cd,
            const int i_center,
            const std::vector<int>& disk_indices)
    {
        for (const int offset : disk_indices)
        {
            const int ijk = i_center + offset;

            const TF u_on_u = u[ijk] + utrans;
            const TF ftau   = -cd * std::abs(u_on_u);

            ut[ijk] += ftau * u_on_u;
        }
    }

    // TODO: implement disk_avg_velocity:
    // - must be called after initializing disk points, but before force computation
    // Sample velocities at turbine indices 
    // take weighted (use weighted when area fraction on disk is used) average of velocities at point

    // TODO: implement compute_disk_forces:
    /* Computes the axial (thrust) and tangential (rotational) forces on disk
    *
    *   Inputs:
    *   turbine index, disk indices (flattened), disk velocities
    * 
    *   Compute axial thrust force:
    *   -> 0.5 * Ct * U^2 * A_{disk-or-cell}/dx
    * 
    *   Compute radial tangential force:
    *   -> 0.5 * Cp * U^2 * (U/(omega*r)) * cos(theta) * A_{disk-or-cell}/dx
    *   -> 0.5 * Cp * U^2 * (U/(omega*r)) * sin(theta) * A_{disk-or-cell}/dx
    * 
    */

    // TODO: implement compute_disk_power
    /* Computes the power at the rotor disk
    *
    *   -> reference lines 730-752 of windfarm_module_v1.f90 in Fabien's code
    */

    // TODO: implement compute_disk_torque
    /* Computes the torque at the rotor disk
    *
    *   -> reference lines 730-752 of windfarm_module_v1.f90 in Fabien's code
    * 
    */
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
        diameter = input.get_item<TF>("dragdisk", "diameter", "", TF(100));
        height   = input.get_item<TF>("dragdisk", "height", "", TF(100));
        cd       = input.get_item<TF>("dragdisk", "cd", "", TF(0.5));
        xc       = input.get_item<TF>("dragdisk", "xc", "", TF(100));
        yc       = input.get_item<TF>("dragdisk", "yc", "", TF(100));
    }
}

template<typename TF>
DragDisk<TF>::~DragDisk()
{
}

template <typename TF>
void DragDisk<TF>::init()
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

    // TODO: Add loop through the number of turbines, then call find_disk_points at each turbine
    // location and assign each set of turbine location indices to a specific turbine index

    find_disk_points(grid, height, diameter, xc, yc, i_center, disk_indices);
}

#ifndef USECUDA
template<typename TF>
void DragDisk<TF>::exec(Stats<TF>& stats)
{
    if (!sw_disk)
        return;

    if (disk_indices.empty())
        return;

    auto& gd = grid.get_grid_data();

    // TODO: Call actuator disk function that gets mean velocity at the disk itself, use to calculate the velocity used in thrust force
    // Include the wake rotation (tangential force) in the actuator disk

    dragdisk_u<TF>(
        fields.mt.at("u")->fld.data(),  // ut: u-tendency array
        fields.mp.at("u")->fld.data(),  // u:  prognostic u field
        gd.utrans,
        cd,
        i_center,    
        disk_indices);

    stats.calc_tend(*fields.mt.at("u"), tend_name);
}
#endif


#ifdef FLOAT_SINGLE
template class DragDisk<float>;
#else
template class DragDisk<double>;
#endif
