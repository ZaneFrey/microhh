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

#include "dragdisk.h"
#include "master.h"
#include "grid.h"
#include "fields.h"
#include "input.h"
#include "stats.h"
#ifdef USECUDA
#include "cuda_buffer.h"
template<typename TF>
void dragdisk_apply_cuda(cuda_vector<TF>& u, const cuda_vector<int>& indices, TF cd);
#endif
#include <algorithm>
#include <cmath>

namespace
{
template<typename TF>
void find_disk_points(const Grid<TF>& grid, TF height, int sizex, int sizey,
                      int& k_center, int& istart, int& iend, int& jstart, int& jend)
{
    const auto& gd = grid.get_grid_data();

    // Find vertical index closest to desired height
    k_center = gd.kstart;
    TF minz = std::abs(gd.z[gd.kstart] - height);
    for (int k = gd.kstart + 1; k < gd.kend; ++k)
    {
        TF diff = std::abs(gd.z[k] - height);
        if (diff < minz)
        {
            minz = diff;
            k_center = k;
        }
    }

    // Determine disk bounds centred in domain
    int ic = gd.istart + (gd.iend - gd.istart) / 2;
    int jc = gd.jstart + (gd.jend - gd.jstart) / 2;

    istart = std::max(gd.istart, ic - sizex / 2);
    iend   = std::min(gd.iend,   istart + sizex);
    jstart = std::max(gd.jstart, jc - sizey / 2);
    jend   = std::min(gd.jend,   jstart + sizey);
}
} // unnamed namespace

template<typename TF>
DragDisk<TF>::DragDisk(Master& masterin, Grid<TF>& gridin, Fields<TF>& fieldsin, Input& input) :
    master(masterin), grid(gridin), fields(fieldsin)
{
    sizex  = input.get_item<int>("dragdisk", "sizex", "", 0);
    sizey  = input.get_item<int>("dragdisk", "sizey", "", 0);
    height = input.get_item<TF>("dragdisk", "height", "", TF(0));
    cd     = input.get_item<TF>("dragdisk", "cd", "", TF(0.2));
    enabled = sizex > 0 && sizey > 0;
    k_center = 0;
}

template<typename TF>
DragDisk<TF>::~DragDisk()
{
}

template<typename TF>
void DragDisk<TF>::create()
{
    if (!enabled)
        return;

    auto& gd = grid.get_grid_data();

    int istart, iend, jstart, jend;
    find_disk_points(grid, height, sizex, sizey, k_center, istart, iend, jstart, jend);

    const int jj = gd.jstride;
    const int kk = gd.kstride;

    indices.clear();
    for (int j = jstart; j < jend; ++j)
        for (int i = istart; i < iend; ++i)
            indices.push_back(i + j * jj + k_center * kk);
}

template<typename TF>
void DragDisk<TF>::exec(Stats<TF>&, double)
{
    if (!enabled)
        return;
#ifdef USECUDA
    dragdisk_apply_cuda(fields.mp.at("u")->fld_g, indices_g, cd);
#else
    auto& u = fields.mp.at("u")->fld;
    for (int idx : indices)
        u[idx] -= cd * u[idx];
#endif
}

template<typename TF>
void DragDisk<TF>::prepare_device()
{
#ifdef USECUDA
    if (!enabled)
        return;
    indices_g.allocate(indices.size());
    cuda_copy(indices, indices_g);
#endif
}

template<typename TF>
void DragDisk<TF>::clear_device()
{
#ifdef USECUDA
    indices_g.free();
#endif
}

#ifdef FLOAT_SINGLE
template class DragDisk<float>;
#else
template class DragDisk<double>;
#endif
