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

#include "turbine.h"
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
    template<typename TF>
    void find_disk_points(
        const Grid<TF>& grid,
        int turbine_idx,
        TF height,
        TF diameter,
        TF xc,
        TF yc,
        int& i_center,
        std::vector<int>& disk_indices,
        std::vector<TF>& radial_dist,
        std::vector<TF>& angle)
    {
        const auto& gd = grid.get_grid_data();

        // Ensure the disk lies inside the physical domain.
        const TF radius = TF(0.5) * diameter;
        if (yc - radius < TF(0) || yc + radius > gd.ysize ||
            height - radius < TF(0) || height + radius > gd.zsize)
        {
            throw std::runtime_error("Turbine: disk boundaries lie outside the global domain");
        }

        // y-direction ownership check: if the disk (including its radius)
        // does not intersect this rank's y-interval, skip it on this rank.
        const TF dy = gd.dy;
        const TF y_min_center = gd.y[gd.jstart];
        const TF y_max_center = gd.y[gd.jend - 1];
        const TF y_min_edge   = y_min_center - TF(0.5) * dy;
        const TF y_max_edge   = y_max_center + TF(0.5) * dy;

        if (yc + radius < y_min_edge || yc - radius > y_max_edge)
        {
            disk_indices.clear();
            return;
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
        const int j_lo = std::max(j_min, gd.jstart);
        const int j_hi = std::min(j_max, gd.jend - 1);
        const int k_lo = std::max(k_min, gd.kstart);
        const int k_hi = std::min(k_max, gd.kend - 1);

        disk_indices.clear();

        if (j_lo > j_hi || k_lo > k_hi)
            return; // Disk does not intersect this rank in y-z.

        // Extract disk indices on y-z plane that fall entirely within the diameter.
        // Indices are flattened as j * jstride + k * kstride
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
            const TF dz  = gd.z[k] - height;
            const TF dz2 = dz * dz;

            for (int j = j_lo; j <= j_hi; ++j)
            {
                const TF dy_off = gd.y[j] - yc;
                const TF dist2  = dy_off * dy_off + dz2;

                if (dist2 <= eff_radius2)
                {
                    const int index = j * jstride + k * kstride;
                    disk_indices.push_back(index);
                }
            }
        }

        // Calculate radial distance and azimuthal angle for use in disk forcing later
        for (const int offset : disk_indices)
        {
            const int ijk = i_center + offset;
            const int k_idx = offset / kstride;
            const int j_idx = (offset - k_idx * kstride) / jstride;
            const TF offset_y_cells = (gd.y[j_idx] - yc) * gd.dyi;
            const TF offset_z_cells = (gd.z[k_idx] - height) / dz_center;
            const int offset_y = static_cast<int>(std::round(offset_y_cells));
            const int offset_z = static_cast<int>(std::round(offset_z_cells));
            
            // Radial distance of every point in disk region from the hub center
            const TF radial_dist[ijk] = std::sqrt(std::pow(offset_y * gd.dy,2) + std::pow(offset_z * gd.dz,2));

            // Azimuthal angle on tangential direction of every point from hub center
            const TF angle[ijk] = std::atan2(offset_z * gd.dz, offset_y * gd.dy) + TF(M_PI_2);
        }
    
    }

    template<typename TF>
    void calc_disk_velocity(
        const int turbine_idx,
        const TF* const restrict u,
        const TF utrans,
        const TF diameter,
        const TF tsr,
        const int i_center
        const std::vector<int>& disk_indices,
        const std::vector<TF>& disk_vel,
        const std::vector<TF>& omega)
    {
        if (disk_indices.empty())
        {
            disk_vel[turbine_idx] = TF(0);
            omega[turbine_idx]    = TF(0);
            return;
        }
        
        // Calculate average rotor velocity at rotor disk
        TF sum = TF(0);
        for (const int offset : disk_indices)
        {
            const int ijk = i_center + offset;
            sum += u[ijk] + utrans;
        }

        disk_vel[turbine_idx] = sum / static_cast<TF>(disk_indices.size());     // disk-average streamwise velocity
        omega[turbine_idx] = tsr * disk_vel[turbine_idx] / (diameter / 2);      // angular velocity
    }

   template<typename TF>
   void calc_disk_forces(
        const int turbine_idx,
        TF* const restrict ut, TF* const restrict vt, TF* const restrict wt,
        const TF* const restrict u, const TF utrans,
        const TF cp, const TF ct,
        const TF area_cell,
        const std::vector<int>& disk_indices,
        const std::vector<TF>& disk_vel, const std::vector<TF>& omega,
        const std::vector<TF>& angle, const std::vector<TF>& radial_dist,
        const int i_center, const int jstride, const int kstride)
    {
        const int ii = 1;
        const int jj = jstride;
        const int kk = kstride;

        for (const int offset : disk_indices)
        {
            const int ijk = i_center + offset;

            // Axial (thrust) forcing
            const TF ftx = -TF(0.5) * ct * (disk_vel[turbine_idx] * disk_vel[turbine_idx]) * area_cell;

            ut[ijk] += ftx;

            // Wake swirl (tangential) forcing projected onto staggered v and w grids.
            if (radial_dist[ijk] != TF(0))
            {
                const TF fty = -TF(0.5) * cp * (disk_vel[turbine_idx] * disk_vel[turbine_idx]) * (u[ijk] + utrans) / 
                                (omega[turbine_idx] * radial_dist[ijk]) * std::cos(angle[ijk]) * area_cell;
                const TF ftz = -TF(0.5) * cp * (disk_vel[turbine_idx] * disk_vel[turbine_idx]) * (u[ijk] + utrans) / 
                                (omega[turbine_idx] * radial_dist[ijk]) * std::sin(angle[ijk]) * area_cell;

                // Spread tangential forces to v, w tendencies
                const TF fty_share = fty * TF(0.25);
                vt[ijk]           += fty_share;
                vt[ijk - ii]      += fty_share;
                vt[ijk - ii + jj] += fty_share;
                vt[ijk + jj]      += fty_share;

                const TF ftz_share = ftz * TF(0.25);
                wt[ijk]           += ftz_share;
                wt[ijk - ii]      += ftz_share;
                wt[ijk - ii + kk] += ftz_share;
                wt[ijk + kk]      += ftz_share;
            }    
        }
    }

    template<typename TF>
    void calc_disk_power_torque(
        const int turbine_idx,
        const TF* const restrict u, const TF utrans,
        const TF cp,
        const TF area_cell,
        const int i_center,
        const std::vector<int>& disk_indices,
        const std::vector<TF>& disk_vel, const std::vector<TF>& omega,
        std::vector<TF>& power, std::vector<TF>& torque)
    {
        power[turbine_idx]  = TF(0);
        torque[turbine_idx] = TF(0);

        for (const int offset : disk_indices)
        {
            const int ijk = i_center + offset;

            // Calculate the turbine power
            const TF pt = TF(0.5) * cp * (u[ijk] + utrans) *
                        (disk_vel[turbine_idx] * disk_vel[turbine_idx]) * area_cell;
            power[turbine_idx] = power[turbine_idx] + pt;
            
            // Calculate the turbine torque
                const TF tt = pt / omega[turbine_idx];
                torque[turbine_idx] = torque[turbine_idx] + tt;
        }
    }

    // TODO: implement tower_nacelle_drag:
    /* Represents turbine tower and nacelle geometry as subgrid canopy-style elements with drag forcing
    *
    *   ->  Towers are represented with TAD (tower area density, analogous to PAD),
    *       with units of m^2/m^3 (frontal area/cell volume)
    *   ->  Tower diameter tapers with height, will result in decreasing TAD
    *   ->  Uses Cd from tabulated data of flow across cylinders, Cd~0.5 for Re_D~1e6
    *   ->  Takes turbine disk center from find_disk_points
    */
    template<typename TF>
    void tower_nacelle_drag(
            TF* const restrict ut,
            TF* const restrict vt,
            TF* const restrict vw,
            const TF* const restrict u,
            const TF* const restrict v,
            const TF* const restrict w,
            const TF* const restrict tad,
            const TF utrans,
            const TF vtrans,
            const TF tower_diameter,
            const TF height,
            const TF cd,
            const int i_center,
            const int j_center,
            const int k_center)
    {}
} 

template<typename TF>
Turbine<TF>::Turbine(
        Master& masterin, Grid<TF>& gridin, Fields<TF>& fieldsin, Input& input) :
        master(masterin), grid(gridin), fields(fieldsin), 
        field3d_operators(masterin, gridin, fieldsin), field3d_io(masterin, gridin)
{
    sw_adm      = input.get_item<bool>("turbine", "swadm", "", false);
    sw_windfarm = input.get_item<bool>("turbine", "swwindfarm", "", false);

    if (sw_adm)
    {
        if (sw_windfarm)
        {
            diameter = input.get_item<TF>("turbine", "diameter", "", TF(100));
            height   = input.get_item<TF>("turbine", "height",   "", TF(100));
            cp       = input.get_item<TF>("turbine", "cp",       "", TF(0.45));
            ct       = input.get_item<TF>("turbine", "ct",       "", TF(1.333));
            tsr      = input.get_item<TF>("turbine", "tsr",      "", TF(8));
            // TODO: add read from list to get xc and yc locations from list for multiple turbines
        }
        else
        {
            turbine_idx = 1;
            diameter = input.get_item<TF>("turbine", "diameter", "", TF(100));
            height   = input.get_item<TF>("turbine", "height",   "", TF(100));
            cp       = input.get_item<TF>("turbine", "cp",       "", TF(0.45));
            ct       = input.get_item<TF>("turbine", "ct",       "", TF(1.333));
            tsr      = input.get_item<TF>("turbine", "tsr",      "", TF(8));
            xc       = input.get_item<TF>("turbine", "xc",       "", TF(100));
            yc       = input.get_item<TF>("turbine", "yc",       "", TF(100));
        }
    }

}

template<typename TF>
Turbine<TF>::~Turbine()
{
}

template <typename TF>
void Turbine<TF>::init()
{
    if (!sw_adm)
        return;

    auto& gd = grid.get_grid_data();

}

template<typename TF>
void Turbine<TF>::create(
        Input& inputin, Stats<TF>& stats)
{
    if (!sw_adm)
        return;

    // TODO: Add loop through the number of turbines, then call find_disk_points at each turbine
    // location and assign each set of turbine location indices to a specific turbine index

    find_disk_points(
        grid,
        turbine_idx,
        height,
        diameter,
        xc,
        yc,
        i_center,
        disk_indices,
        radial_dist,
        angle);
}

#ifndef USECUDA
template<typename TF>
void Turbine<TF>::exec(Stats<TF>& stats)
{
    if (!sw_adm)
        return;

    if (disk_indices.empty())
        return;

    auto& gd = grid.get_grid_data();

    // Calculate disk area and equivalent cell areas
    const TF area_disk          = TF(M_PI) * diameter * diameter / TF(4);
    const std::size_t disk_npts = disk_indices.size();
    const TF area_cell          = area_disk / static_cast<TF>(disk_npts);

    // Compute average disk velocity
    TF disk_vel[turbine_idx];
    calc_disk_velocity<TF>(
        turbine_idx,
        fields.mp.at("u")->fld.data(),
        gd.utrans,
        diameter,
        tsr,
        i_center,
        disk_indices,
        disk_vel,
        omega);

    // Compute forces at each disk (axial and tangential)
    calc_disk_forces<TF>(
        turbine_idx,
        fields.mp.at("u")->fld.data(),
        fields.mp.at("v")->fld.data(),
        fields.mp.at("w")->fld.data(),
        fields.mp.at("u")->fld.data(),
        gd.utrans,
        cp,
        ct,
        area_cell,
        disk_indices,
        disk_vel[turbine_idx], omega[turbine_idx],
        angle, radial_dist,
        i_center, gd.icells, gd.ijcells);

    // Disk power, integrated over the disk.
    calc_disk_power_torque<TF>(
        turbine_idx,
        fields.mp.at("u")->fld.data(), gd.utrans,
        cp,
        area_cell,
        i_center,
        disk_indices,
        disk_vel[turbine_idx], omega[turbine_idx],
        power[turbine_idx], torque[turbine_idx]);
}
#endif


#ifdef FLOAT_SINGLE
template class Turbine<float>;
#else
template class Turbine<double>;
#endif
