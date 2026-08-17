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
#include <stdexcept>

#include "turbine.h"
#include "grid.h"
#include "fields.h"
#include "stats.h"
#include "tools.h"

#include "turbine_kernels.cuh"
#include "cuda_launcher.h"
#include "cuda_tiling.h"

namespace
{
    
}

#ifdef USECUDA
template<typename TF>
void Turbine<TF>::prepare_device()
{
    // No persistent device-side turbine data is required at the moment.
}

template<typename TF>
void Turbine<TF>::clear_device()
{
    // No persistent device-side turbine data has been allocated.
}

template<typename TF>
void Turbine<TF>::exec(const double dt, Stats<TF>& stats)
{
    (void)stats;

    if (!sw_adm)
        return;

    // If no turbines were created on this rank, there is nothing to do.
    if (disk_indices.empty())
        return;

    if (sw_tower)
    {
        throw std::runtime_error("Turbine: swtower is enabled, but tower/nacelle drag is not implemented for CUDA builds");
    }

    const auto& gd = grid.get_grid_data();

    // Weighting function for time-filtering
    const TF dt_TF = static_cast<TF>(dt);
    if (sw_timefilter)
    {
        const TF ratio = dt_TF / disk_mem_time;
        disk_weight         = ratio / (TF(1) + ratio);
    }
    else
    {
        // Disable time-filtering: use instantaneous disk-space average.
        disk_weight = TF(1);
    }

    // Geometrical disk area (same as CPU implementation).
    const TF area_disk = TF(M_PI) * diameter * diameter / TF(4);

    for (int it = 0; it < nturbine; ++it)
    {
        // CPU create() leaves disk_indices[it] empty if this turbine does
        // not intersect the local sub-domain in y-z.
        if (disk_indices[it].empty())
            continue;

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

        // Compute effective radius following the CPU implementation
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

        // Restrict the launch to the rotor i-plane at i_center[it].
        Grid_layout grid_layout_disk = {
                i_center[it], i_center[it] + 1,
                gd.jstart,    gd.jend,
                gd.kstart,    gd.kend,
                gd.istride,
                gd.jstride,
                gd.kstride};

        // --------------------------------------------------------------------
        // 1. Compute disk-averaged velocity 
        // --------------------------------------------------------------------
        TF sum_vel_host    = TF(0);
        TF sum_weight_host = TF(0);

        TF* sum_vel_dev    = nullptr;
        TF* sum_weight_dev = nullptr;

        cuda_safe_call(cudaMalloc(&sum_vel_dev,    sizeof(TF)));
        cuda_safe_call(cudaMalloc(&sum_weight_dev, sizeof(TF)));
        cuda_safe_call(cudaMemcpy(sum_vel_dev,    &sum_vel_host,    sizeof(TF), cudaMemcpyHostToDevice));
        cuda_safe_call(cudaMemcpy(sum_weight_dev, &sum_weight_host, sizeof(TF), cudaMemcpyHostToDevice));

        launch_grid_kernel<Turbine_kernels::calc_disk_velocity_g<TF>>(
                grid_layout_disk,
                fields.mp.at("u")->fld_g,
                gd.utrans,
                gd.y_g,
                gd.z_g,
                yc[it],
                height,
                eff_radius2,
                sum_vel_dev,
                sum_weight_dev);
        cuda_check_error();
        cudaDeviceSynchronize();

        cuda_safe_call(cudaMemcpy(&sum_vel_host,    sum_vel_dev,    sizeof(TF), cudaMemcpyDeviceToHost));
        cuda_safe_call(cudaMemcpy(&sum_weight_host, sum_weight_dev, sizeof(TF), cudaMemcpyDeviceToHost));

        cudaFree(sum_vel_dev);
        cudaFree(sum_weight_dev);

        if (sum_weight_host > TF(0))
        {
            const TF disk_space_avg = sum_vel_host / sum_weight_host;

            // Apply optional time-filtering to the disk-averaged velocity.
            disk_vel[it] = (TF(1) - disk_weight) * disk_vel[it] + disk_weight * disk_space_avg;
            omega[it]    = tsr * disk_vel[it] / (diameter / TF(2));
        }
        else
        {
            disk_vel[it] = TF(0);
            omega[it]    = TF(0);
            continue; // this turbine does not intersect this rank in practice
        }

        // --------------------------------------------------------------------
        // 2. Apply axial and tangential forces on the rotor disk
        // --------------------------------------------------------------------
        launch_grid_kernel<Turbine_kernels::calc_disk_forces_g<TF>>(
                grid_layout_disk,
                fields.mt.at("u")->fld_g.view(),
                fields.mt.at("v")->fld_g.view(),
                fields.mt.at("w")->fld_g.view(),
                fields.mp.at("u")->fld_g,
                gd.utrans,
                gd.y_g,
                gd.z_g,
                yc[it],
                height,
                eff_radius2,
                cp,
                ct,
                gd.dx,
                disk_vel[it],
                omega[it]);
        cuda_check_error();
        cudaDeviceSynchronize();

        // --------------------------------------------------------------------
        // 3. Compute turbine power and torque on the host side
        // --------------------------------------------------------------------
        if (disk_vel[it] > TF(0) && omega[it] > TF(0))
        {
            power[it]  = TF(0.5) * cp * disk_vel[it] * disk_vel[it] * disk_vel[it] * area_disk;
            torque[it] = power[it] / omega[it];
        }
        else
        {
            power[it]  = TF(0);
            torque[it] = TF(0);
        }
    }
}
#endif


#ifdef FLOAT_SINGLE
template class Turbine<float>;
#else
template class Turbine<double>;
#endif
