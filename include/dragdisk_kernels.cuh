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

#ifndef DRAGDISK_KERNELS_CUH
#define DRAGDISK_KERNELS_CUH

#include <cmath>

#include "cuda_tiling.h"

namespace DragDisk_kernels
{
    template<typename TF>
    struct dragdisk_u_g
    {
        DEFINE_GRID_KERNEL("dragdisk::dragdisk_u_g", 0)

        template <typename Level>
        CUDA_DEVICE
        void operator()(
                Grid_layout g, const int i, const int j, const int k, const Level level,
                TF* const __restrict__ ut,
                const TF* const __restrict__ u,
                const TF utrans,
                const TF cd,
                const int i_center,
                const int j_center,
                const TF radius_cells_sq)
        {
            (void)level;

            const TF di = TF(i - i_center);
            const TF dj = TF(j - j_center);
            if (di * di + dj * dj > radius_cells_sq)
                return;

            const int ijk = g(i, j, k);
            const TF u_on_u = u[ijk] + utrans;
            const TF ftau = -cd * fabs(u_on_u);

            ut[ijk] += ftau * u_on_u;
        }
    };
}

#endif // DRAGDISK_KERNELS_CUH
