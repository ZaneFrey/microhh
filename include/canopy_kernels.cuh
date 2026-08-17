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

#ifndef CANOPY_KERNELS_CUH
#define CANOPY_KERNELS_CUH

#include <cmath>

#include "cuda_tiling.h"

namespace Canopy_kernels
{
    template<typename TF>
    struct canopy_drag_u_g
    {
        DEFINE_GRID_KERNEL("canopy::canopy_drag_u_g", 0)

        template <typename Level>
        CUDA_DEVICE
        void operator()(
                Grid_layout g, const int i, const int j, const int k, const Level level,
                TF* const __restrict__ ut,
                const TF* const __restrict__ u,
                const TF* const __restrict__ v,
                const TF* const __restrict__ w,
                const TF* const __restrict__ pad,
                const TF utrans,
                const TF vtrans,
                const TF cd)
        {
            (void)level;

            const int ii = 1;
            const int jj = g.jstride;
            const int kk = g.kstride;

            const int ijk = g(i, j, k);

            const TF u_on_u = u[ijk] + utrans;
            const TF v_on_u = TF(0.25) * (v[ijk] + v[ijk-ii] + v[ijk-ii+jj] + v[ijk+jj]) + vtrans;
            const TF w_on_u = TF(0.25) * (w[ijk] + w[ijk-ii] + w[ijk-ii+kk] + w[ijk+kk]);

            const TF speed = sqrt(u_on_u*u_on_u + v_on_u*v_on_u + w_on_u*w_on_u);
            const TF ftau = -cd * pad[k] * speed;

            ut[ijk] += ftau * u_on_u;
        }
    };

    template<typename TF>
    struct canopy_drag_v_g
    {
        DEFINE_GRID_KERNEL("canopy::canopy_drag_v_g", 0)

        template <typename Level>
        CUDA_DEVICE
        void operator()(
                Grid_layout g, const int i, const int j, const int k, const Level level,
                TF* const __restrict__ vt,
                const TF* const __restrict__ u,
                const TF* const __restrict__ v,
                const TF* const __restrict__ w,
                const TF* const __restrict__ pad,
                const TF utrans,
                const TF vtrans,
                const TF cd)
        {
            (void)level;

            const int ii = 1;
            const int jj = g.jstride;
            const int kk = g.kstride;

            const int ijk = g(i, j, k);

            const TF u_on_v = TF(0.25) * (u[ijk] + u[ijk+ii] + u[ijk-jj] + u[ijk+ii-jj]) + utrans;
            const TF v_on_v = v[ijk] + vtrans;
            const TF w_on_v = TF(0.25) * (w[ijk] + w[ijk+kk] + w[ijk-jj] + w[ijk+kk-jj]);

            const TF speed = sqrt(u_on_v*u_on_v + v_on_v*v_on_v + w_on_v*w_on_v);
            const TF ftau = -cd * pad[k] * speed;

            vt[ijk] += ftau * v_on_v;
        }
    };

    template<typename TF>
    struct canopy_drag_w_g
    {
        DEFINE_GRID_KERNEL("canopy::canopy_drag_w_g", 0)

        template <typename Level>
        CUDA_DEVICE
        void operator()(
                Grid_layout g, const int i, const int j, const int k, const Level level,
                TF* const __restrict__ wt,
                const TF* const __restrict__ u,
                const TF* const __restrict__ v,
                const TF* const __restrict__ w,
                const TF* const __restrict__ padh,
                const TF utrans,
                const TF vtrans,
                const TF cd)
        {
            (void)level;

            const int ii = 1;
            const int jj = g.jstride;
            const int kk = g.kstride;

            const int ijk = g(i, j, k);

            const TF u_on_w = TF(0.25) * (u[ijk] + u[ijk+ii] + u[ijk-kk] + u[ijk+ii-kk]) + utrans;
            const TF v_on_w = TF(0.25) * (v[ijk] + v[ijk+jj] + v[ijk-kk] + v[ijk+jj-kk]) + vtrans;
            const TF w_on_w = w[ijk];

            const TF speed = sqrt(u_on_w*u_on_w + v_on_w*v_on_w + w_on_w*w_on_w);
            const TF ftau = -cd * padh[k] * speed;

            wt[ijk] += ftau * w_on_w;
        }
    };
}

#endif // CANOPY_KERNELS_CUH
