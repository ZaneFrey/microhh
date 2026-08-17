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

#include "average.h"
#include "cuda_buffer.h"
#include "dump.h"
#include "grid.h"
#include "master.h"
#include "tools.h"

namespace
{
    template<typename TF> __global__
    void accumulate_kernel_g(
            TF* const average, const TF* const field, const TF dt, const int ncells)
    {
        const int n = blockIdx.x*blockDim.x + threadIdx.x;

        if (n < ncells)
            average[n] += field[n] * dt;
    }

    template<typename TF> __global__
    void scale_g(TF* const field, const TF fac, const int ncells)
    {
        const int n = blockIdx.x*blockDim.x + threadIdx.x;

        if (n < ncells)
            field[n] *= fac;
    }

    template<typename TF> __global__
    void set_to_zero_g(TF* const field, const int ncells)
    {
        const int n = blockIdx.x*blockDim.x + threadIdx.x;

        if (n < ncells)
            field[n] = TF(0);
    }
}

template<typename TF>
void Average<TF>::prepare_device()
{
    auto& gd = grid.get_grid_data();

    for (auto& it : average_fields)
        it.second.fld_g.resize(gd.ncells);

    forward_device();
}

template<typename TF>
void Average<TF>::clear_device()
{
    for (auto& it : average_fields)
        it.second.fld_g.free();
}

template<typename TF>
void Average<TF>::forward_device()
{
    auto& gd = grid.get_grid_data();

    for (auto& it : average_fields)
        cuda_copy(it.second.fld.data(), it.second.fld_g.data(), gd.ncells);
}

template<typename TF>
void Average<TF>::accumulate_g(const std::string& name, const TF* const data, const TF dt)
{
    auto& gd = grid.get_grid_data();
    auto& average_field = average_fields.at(name);

    const int blocki = 256;
    const int gridi = gd.ncells/blocki + (gd.ncells%blocki > 0);

    accumulate_kernel_g<TF><<<gridi, blocki>>>(average_field.fld_g, data, dt, gd.ncells);
    cuda_check_error();
}

template<typename TF>
void Average<TF>::reset_fields()
{
    auto& gd = grid.get_grid_data();
    const int blocki = 256;
    const int gridi = gd.ncells/blocki + (gd.ncells%blocki > 0);

    for (auto& it : average_fields)
    {
        set_to_zero_g<TF><<<gridi, blocki>>>(it.second.fld_g, gd.ncells);
        cuda_check_error();
        std::fill(it.second.fld.begin(), it.second.fld.end(), TF(0));
    }
}

template<typename TF>
void Average<TF>::write_field(const std::string& name, int iotime)
{
    auto& gd = grid.get_grid_data();
    auto& average_field = average_fields.at(name);
    const int blocki = 256;
    const int gridi = gd.ncells/blocki + (gd.ncells%blocki > 0);
    const TF fac = TF(1) / TF(average_time);

    scale_g<TF><<<gridi, blocki>>>(average_field.fld_g, fac, gd.ncells);
    cuda_check_error();

    cuda_copy(average_field.fld_g.data(), average_field.fld.data(), gd.ncells);
    dump.save_data(average_field.fld.data(), name + "_avg", iotime);
}

template<typename TF>
void Average<TF>::save(const int iotime)
{
    if (!swaverage || !has_fields())
        return;

    for (auto& it : average_fields)
    {
        auto& gd = grid.get_grid_data();
        cuda_copy(it.second.fld_g.data(), it.second.fld.data(), gd.ncells);
        dump.save_data(it.second.fld.data(), it.first + "_avg_acc", iotime);
    }

    save_average_time(iotime);
}

#ifdef FLOAT_SINGLE
template class Average<float>;
#else
template class Average<double>;
#endif
