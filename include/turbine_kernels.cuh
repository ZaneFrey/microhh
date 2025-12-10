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

#ifndef TURBINE_KERNELS_CUH
#define TURBINE_KERNELS_CUH

#include <cmath>

#include "cuda_tiling.h"

namespace Turbine_kernels
{
    template<typename TF>
    struct calc_disk_velocity_g
    {
        DEFINE_GRID_KERNEL("turbine::calc_disk_velocity_g", 0)

        template <typename Level>
        CUDA_DEVICE
        void operator()(
                Grid_layout g, const int i, const int j, const int k, const Level level,
                const TF* const __restrict__ u,
                const TF utrans,
                const TF* const __restrict__ y,
                const TF* const __restrict__ z,
                const TF yc,
                const TF height,
                const TF eff_radius2,
                TF* const __restrict__ sum_vel,
                TF* const __restrict__ sum_weight)
        {
            (void)level;

            const TF dy    = y[j] - yc;
            const TF dz    = z[k] - height;
            const TF dist2 = dy * dy + dz * dz;

            if (dist2 > eff_radius2)
                return;

            const int ijk    = g(i, j, k);
            const TF u_local = u[ijk] + utrans;

            atomicAdd(sum_vel,    u_local);
            atomicAdd(sum_weight, TF(1));
        }
    };

    template<typename TF>
    struct calc_disk_forces_g
    {
        DEFINE_GRID_KERNEL("turbine::calc_disk_forces_g", 0)

        template <typename Level>
        CUDA_DEVICE
        void operator()(
                Grid_layout g, const int i, const int j, const int k, const Level level,
                TF* const __restrict__ ut,
                TF* const __restrict__ vt,
                TF* const __restrict__ wt,
                const TF* const __restrict__ u,
                const TF utrans,
                const TF* const __restrict__ y,
                const TF* const __restrict__ z,
                const TF yc,
                const TF height,
                const TF eff_radius2,
                const TF cp,
                const TF ct,
                const TF dx,
                const TF disk_vel,
                const TF omega)
        {
            (void)level;

            const TF dy    = y[j] - yc;
            const TF dz    = z[k] - height;
            const TF dist2 = dy * dy + dz * dz;

            if (dist2 > eff_radius2)
                return;

            if (disk_vel == TF(0) || omega == TF(0))
                return;

            const int ijk = g(i, j, k);

            const int ii = 1;
            const int jj = g.jstride;
            const int kk = g.kstride;

            // Axial (thrust) forcing (constant over the disk for a given turbine).
            const TF ftx = -TF(0.5) * ct * (disk_vel * disk_vel) / dx;
            ut[ijk] += ftx;

            const TF r2 = dist2;
            if (r2 == TF(0))
                return;

            const TF r   = std::sqrt(r2);
            const TF phi = std::atan2(dz, dy) + TF(M_PI_2);

            const TF u_local = u[ijk] + utrans;
            const TF factor  = -TF(0.5) * cp * (disk_vel * disk_vel) * u_local /
                               (omega * r) / dx;

            const TF fty = factor * std::cos(phi);
            const TF ftz = factor * std::sin(phi);

            // Spread tangential forces to v, w tendencies using the same 4-point stencil
            // as in the CPU implementation.
            const TF fty_share = fty * TF(0.25);
            atomicAdd(&vt[ijk],           fty_share);
            atomicAdd(&vt[ijk - ii],      fty_share);
            atomicAdd(&vt[ijk - ii + jj], fty_share);
            atomicAdd(&vt[ijk + jj],      fty_share);

            const TF ftz_share = ftz * TF(0.25);
            atomicAdd(&wt[ijk],           ftz_share);
            atomicAdd(&wt[ijk - ii],      ftz_share);
            atomicAdd(&wt[ijk - ii + kk], ftz_share);
            atomicAdd(&wt[ijk + kk],      ftz_share);
        }
    };

    template<typename TF>
    struct calc_disk_power_torque_g
    {
        DEFINE_GRID_KERNEL("turbine::calc_disk_power_torque_g", 0)

        template <typename Level>
        CUDA_DEVICE
        void operator()(
                Grid_layout g, const int i, const int j, const int k, const Level level,
                TF* const __restrict__ st,
                const TF* const __restrict__ s,
                const TF* const __restrict__ wls,
                const TF* const __restrict__ dzhi)
        {
            // do something
        }
    };
}

#endif // TURBINE_KERNELS_CUH
