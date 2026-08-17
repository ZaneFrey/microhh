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
#include "finite_difference.h"
#include "input.h"
#include "master.h"
#include "netcdf_interface.h"
#include "stats.h"

namespace
{
    template<typename TF>
    void find_disk_points(
        const Grid<TF>& grid,
        TF height,
        TF diameter,
        TF xc,
        TF yc,
        int& i_center,
        int& k_hub,
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
        int    k_hub = gd.kstart;
        TF     k_best   = std::abs(gd.z[k_hub] - height);
        for (int k = gd.kstart + 1; k < gd.kend; ++k)
        {
            const TF diff = std::abs(gd.z[k] - height);
            if (diff < k_best)
            {
                k_best   = diff;
                k_hub = k;
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
        if (k_hub > gd.kstart && k_hub < gd.kend - 1)
            dz_center = TF(0.5) * (gd.z[k_hub + 1] - gd.z[k_hub - 1]); // interior cell
        else if (k_hub < gd.kend - 1)
            dz_center = gd.z[k_hub + 1] - gd.z[k_hub]; // bottom edge case
        else
            dz_center = gd.z[k_hub] - gd.z[k_hub - 1]; // top edge case

        int k_extent = static_cast<int>(std::ceil(radius / dz_center));

        const int j_min = j_center - j_extent;
        const int j_max = j_center + j_extent;
        const int k_min = k_hub - k_extent;
        const int k_max = k_hub + k_extent;

        // Intersect the disk bounding box with the local sub-domain.
        const int j_lo = std::max(j_min, gd.jstart);
        const int j_hi = std::min(j_max, gd.jend - 1);
        const int k_lo = std::max(k_min, gd.kstart);
        const int k_hi = std::min(k_max, gd.kend - 1);

        disk_indices.clear();
        radial_dist.clear();
        angle.clear();

        if (j_lo > j_hi || k_lo > k_hi)
            return; // Disk does not intersect this rank in y-z.

        // Extract disk indices on y-z plane that fall entirely within the diameter.
        // Indices are flattened as j * jstride + k * kstride
        const TF cell_diag = std::sqrt(dy * dy + dz_center * dz_center);
        const TF eff_radius = std::max(TF(0), radius - TF(0.5) * cell_diag);
        const TF eff_radius2 = eff_radius * eff_radius;

        const int ny = j_hi - j_lo + 1;
        const int nz = k_hi - k_lo + 1;
        const int nmax = ny * nz;
        disk_indices.reserve(nmax);
        radial_dist.reserve(nmax);
        angle.reserve(nmax);

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

                    // Radial distance and azimuthal angle for each point in the disk region
                    const TF r   = std::sqrt(dist2);
                    const TF phi = std::atan2(dz, dy_off) + TF(M_PI_2);

                    radial_dist.push_back(r);
                    angle.push_back(phi);
                }
            }
        }
    }

    template<typename TF>
    void calc_disk_velocity(
        const TF* const restrict u,
        const TF utrans,
        const TF diameter,
        const TF tsr,
        const TF disk_vel_old,
        const TF disk_weight,
        const int i_center,
        const std::vector<int>& disk_indices,
        TF& disk_vel,
        TF& omega)
    {
        if (disk_indices.empty())
        {
            disk_vel = TF(0);
            omega    = TF(0);
            return;
        }

        TF sum = TF(0);
        for (const int offset : disk_indices)
        {
            const int ijk = i_center + offset;
            sum += u[ijk] + utrans;
        }

        const TF disk_space_avg = sum / static_cast<TF>(disk_indices.size());

        disk_vel = ((TF(1) - disk_weight) * disk_vel_old) + (disk_weight * disk_space_avg);

        omega = tsr * disk_vel / (diameter / TF(2));
    }

   template<typename TF>
   void calc_disk_forces(
        TF* const restrict ut, TF* const restrict vt, TF* const restrict wt,
        const TF* const restrict u, const TF utrans,
        const TF cp, const TF ct,
        const TF dx,
        const std::vector<int>& disk_indices,
        const TF disk_vel, const TF omega,
        const std::vector<TF>& angle,
        const std::vector<TF>& radial_dist,
        const int i_center, const int jstride, const int kstride)
    {
        const int ii = 1;
        const int jj = jstride;
        const int kk = kstride;

        // low omega check
        if (disk_vel == 0)
        {
            TF ftx = TF(0);
            TF fty = TF(0);
            TF ftz = TF(0);
        }
        else
        {
            for (std::size_t p = 0; p < disk_indices.size(); ++p)
            {
                const int offset = disk_indices[p];
                const int ijk = i_center + offset;
                const TF r = radial_dist[p];
                const TF phi = angle[p];

                // Axial (thrust) forcing
                const TF ftx = -TF(0.5) * ct * (disk_vel * disk_vel) / dx;

                ut[ijk] += ftx;

                // Wake swirl (tangential) forcing projected onto staggered v and w grids.
                if (r != TF(0))
                {
                    const TF fty = -TF(0.5) * cp * (disk_vel * disk_vel) * (u[ijk] + utrans) / 
                                    (omega * r) * std::cos(phi) / dx;
                    const TF ftz = -TF(0.5) * cp * (disk_vel * disk_vel) * (u[ijk] + utrans) / 
                                    (omega * r) * std::sin(phi) / dx;

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
    }

    template<typename TF>
    void calc_disk_power_torque(
        const TF* const restrict u, const TF utrans,
        const TF cp,
        const TF area_cell,
        const int i_center,
        const std::vector<int>& disk_indices,
        const TF disk_vel, const TF omega,
        TF& power, TF& torque)
    {
        power  = TF(0);
        torque = TF(0);

        for (const int offset : disk_indices)
        {
            const int ijk = i_center + offset;

            // Calculate the turbine power
            const TF pt = TF(0.5) * cp * (u[ijk] + utrans) * (disk_vel * disk_vel) * area_cell;
            power = power + pt;
            
            // Calculate the turbine torque
                const TF tt = pt / omega;
                torque = torque + tt;
        }
    }

    // TODO: implement tower_nacelle_drag:
    /* Represents turbine tower and nacelle geometry as subgrid canopy-style elements with drag forcing
    *
    *   ->  Towers are represented with TAD (tower area density, analogous to PAD),
    *       with units of m^2/m^3 (frontal area/cell volume)
    *   ->  Tower is assumed to be u constant-diameter cylinder with height, reducing TAD to diameter/(dx*dy)
    *   ->  Uses Cd from tabulated data of flow across cylinders, Cd~0.5 for Re_D~1e6
    *   ->  Checks whether tower diameter is greater than grid spacing dx/dy. If greater than, 
    *       throw error (will develop later)
    *   ->  Inputs the disk physical center xc/yc and hub height index in z
    *   ->  Loops through bottom z (k) to hub height (k_hub)
    *   ->  Computed TAD
    *   ->  Splits ut and vt calculation to appropriate locations on staggered grid
    *   ->  Apply drag forcing to u and v components vertically, centered on (xc,yc)
    *   ->  Input tower_mem_time, analogous to disk_mem_time but for tower time-filtering
    *   ->  Calculate exponential time-filtering, similar to in turbine disk velocity calculation,
    *       (tower_weight should be calculated in Turbine::exec, similar to calculation of disk_weight)
    *   ->  Update ut vt (tendencies) according to tower drag
    */
    // template<typename TF>
    // void tower_nacelle_drag(
    //         TF* const restrict ut,
    //         TF* const restrict vt,
    //         const TF* const restrict u,
    //         const TF* const restrict v,
    //         const TF utrans,
    //         const TF vtrans,
    //         const TF tower_diameter,
    //         const TF cd,
    //         const TF tower_weight,
    //         const Grid_data<TF>& gd,
    //         const TF xc,
    //         const TF yc,
    //         const int k_hub)
    // {
    //     // Calculate tower area density:
    //     TF tad = tower_diameter / (gd.dx * gd.dy);
        
    //     for (int k = 0; k < k_hub; ++k)
    //     {
    //         // Calculate time-filtered tower velocity
    //         tower_vel = ((TF(1) - tower_weight) * tower_vel_old) + (tower_weight * tower_speed);

    //     }
    // }
} 

template<typename TF>
Turbine<TF>::Turbine(
        Master& masterin, Grid<TF>& gridin, Fields<TF>& fieldsin, Input& input) :
        master(masterin), grid(gridin), fields(fieldsin), 
        field3d_operators(masterin, gridin, fieldsin), field3d_io(masterin, gridin)
{
    sw_adm      = input.get_item<bool>("turbine", "swadm", "", false);
    sw_tower    = input.get_item<bool>("turbine", "swtower", "", false);
    sw_dynyaw   = input.get_item<bool>("turbine", "swdynyaw", "", false);

    if (sw_adm)
    {
        diameter = input.get_item<TF>("turbine", "diameter", "", TF(100));
        height   = input.get_item<TF>("turbine", "height",   "", TF(100));
        cp       = input.get_item<TF>("turbine", "cp",       "", TF(0.8));
        ct       = input.get_item<TF>("turbine", "ct",       "", TF(1.333));
        tsr      = input.get_item<TF>("turbine", "tsr",      "", TF(8));
        xc       = input.get_list<TF>("turbine", "xc", "", std::vector<TF>());
        yc       = input.get_list<TF>("turbine", "yc", "", std::vector<TF>());

        if (xc.empty())
            throw std::runtime_error("Turbine: [turbine][xc] list must not be empty");

        if (xc.size() != yc.size())
            throw std::runtime_error("Turbine: [turbine][xc] and [turbine][yc] lists must have the same length");

        nturbine = static_cast<int>(xc.size());
    
        turbine_starttime   = input.get_item<TF>("turbine","admstarttime","",TF(0));
        disk_mem_time       = input.get_item<TF>("turbine","diskmemtime","",TF(600));
        sw_timefilter       = input.get_item<bool>("turbine","swtimefilter","",true);

        // allocate per‑turbine storage
        i_center.resize(nturbine);
        disk_indices.resize(nturbine);
        radial_dist.resize(nturbine);
        angle.resize(nturbine);
        disk_vel.resize(nturbine);
        omega.resize(nturbine);
        power.resize(nturbine);
        torque.resize(nturbine);

        if (sw_dynyaw)
        {
            yaw_offset      = input.get_list<TF>("turbine", "yawoffset", "", std::vector<TF>());
            scan_dist       = input.get_list<TF>("turbine", "scan_dist", "", std::vector<TF>());
            num_radii       = input.get_item<int>("turbine", "numradii", "", int 50);
            num_azimuths    = input.get_item<int>("turbine", "numazimuths", "", int 50);
            kernel_width    = input.get_item<int>("turbine", "kernelwidth", "", int 3);
        }
    }

    if (sw_tower)
    {
        tower_diameter = input.get_item<TF>("turbine", "towerdiam", "", TF(10));
        tower_mem_time = input.get_item<TF>("turbine", "towermemtime", "", TF(600));
        cd             = input.get_item<TF>("turbine", "cd", "", TF(0.5));
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

    if (sw_tower)
    {
        const auto& gd = grid.get_grid_data();
        tower_speed.resize(nturbine);
        for (int it = 0; it < nturbine; ++it)
            tower_speed[it].assign(gd.kcells, TF(0));

        if (tower_diameter > gd.dx || tower_diameter > gd.dy)
            throw std::runtime_error("Turbine: Tower diameter cannot exceed x/y grid spacing.");
    }

    for (int it = 0; it < nturbine; ++it)
    {
        find_disk_points(
            grid,
            height,
            diameter,
            xc[it],
            yc[it],
            i_center[it],
            k_hub[it],
            disk_indices[it],
            radial_dist[it],
            angle[it]);
    }
}

#ifndef USECUDA
template<typename TF>
void Turbine<TF>::exec(const double dt, Stats<TF>& stats)
{
    if (!sw_adm)
        return;

    if (disk_indices.empty())
        return;

    auto& gd = grid.get_grid_data();


    // Weighting function for time-filerting
    const TF dt_TF = static_cast<TF>(dt);
    if (sw_timefilter)
    {
        const TF ratio_disk = dt_TF / disk_mem_time;
        disk_weight              = ratio_disk / (TF(1) + ratio_disk);

        if (sw_tower)
        {
            const TF ratio_tower = dt_TF / tower_mem_time;
            tower_weight         = ratio_tower / (TF(1) + ratio_tower);
        }
    }
    else // Disable time-filtering: use instantaneous disk-space average
    {
        
        disk_weight = TF(1);
        if (sw_tower)
            tower_weight = TF(1);
    }

    // Calculate disk area and equivalent cell areas
    const TF area_disk = TF(M_PI) * diameter * diameter / TF(4);

    for (int it = 0; it < nturbine; ++it)
    {
        const std::size_t disk_npts = disk_indices[it].size();
        if (disk_npts == 0)
            continue; // this turbine doesn’t intersect this rank

        const TF area_cell = area_disk / static_cast<TF>(disk_npts);

        // Compute average disk velocity
        calc_disk_velocity<TF>(
            fields.mp.at("u")->fld.data(),
            gd.utrans,
            diameter,
            tsr,
            disk_vel[it], // old velocity
            disk_weight,
            i_center[it],
            disk_indices[it],
            disk_vel[it],   // new velocity
            omega[it]);

        // Compute forces at each disk (axial and tangential)
        calc_disk_forces<TF>(
            fields.mt.at("u")->fld.data(),
            fields.mt.at("v")->fld.data(),
            fields.mt.at("w")->fld.data(),
            fields.mp.at("u")->fld.data(),
            gd.utrans,
            cp,
            ct,
            gd.dx,
            disk_indices[it],
            disk_vel[it], 
            omega[it],
            angle[it], 
            radial_dist[it],
            i_center[it], 
            gd.icells,
            gd.ijcells);

        // Disk power, integrated over the rotor
        calc_disk_power_torque<TF>(
            fields.mp.at("u")->fld.data(), 
            gd.utrans,
            cp,
            area_cell,
            i_center[it],
            disk_indices[it],
            disk_vel[it], 
            omega[it],
            power[it], 
            torque[it]);

        if (sw_tower)
        {
            tower_nacelle_drag<TF>(
                fields.mt.at("u")->fld.data(),
                fields.mt.at("v")->fld.data(),
                fields.mp.at("u")->fld.data(),
                fields.mp.at("v")->fld.data(),
                gd.utrans,
                gd.vtrans,
                tower_diameter,
                cd,
                tower_weight,
                gd,
                xc[it],
                yc[it],
                height,
                i_center[it],
                k_hub[it],
                tower_speed[it]);
        }
    }
}
#endif


#ifdef FLOAT_SINGLE
template class Turbine<float>;
#else
template class Turbine<double>;
#endif
