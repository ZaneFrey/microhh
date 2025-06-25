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

#include "turbine.h"
#include "master.h"
#include "grid.h"
#include "fields.h"
#include "input.h"
#include "stats.h"
#include <numeric>
#include <cmath>

template<typename TF>
Turbine<TF>::Turbine(Master& masterin, Grid<TF>& gridin, Fields<TF>& fieldsin,
                     Input& input, TF xin, TF yin, TF hhub_in) :
    master(masterin), grid(gridin), fields(fieldsin)
{
    // Retrieve turbine parameters from input file
    diam          = input.get_item<TF>("turbine", "diam", "");
    hhub          = input.get_item<TF>("turbine", "hhub", "");
    ct            = input.get_item<TF>("turbine", "ct", "");
    cp            = input.get_item<TF>("turbine", "cp", "");
    tsr           = input.get_item<TF>("turbine", "tsr", "");
    swdynyaw      = input.get_item<bool>("turbine", "swdynyaw", "", false);
    yawperiod     = input.get_item<TF>("turbine", "yawperiod", "", TF(0));
    turbstarttime = input.get_item<TF>("turbine", "turbstarttime", "", TF(0));
    swturbstats   = input.get_item<bool>("turbine", "swturbstats", "", false);
    turbstatperiod= input.get_item<TF>("turbine", "turbstatperiod", "", TF(0));

    // Allow override of hub height from layout file
    if (hhub_in > TF(0))
        hhub = hhub_in;

    xpos = xin;
    ypos = yin;
    yaw  = 0;
    next_yaw = turbstarttime;
    // Precompute rotor disk area
    area = M_PI * diam * diam * TF(0.25);
    power = 0;
}

template<typename TF>
Turbine<TF>::~Turbine()
{
}

template<typename TF>
void Turbine<TF>::create()
{
    auto& gd = grid.get_grid_data();

    // Determine vertical index closest to hub height
    k_hub = gd.kstart;
    TF minz = std::abs(gd.z[gd.kstart] - hhub);
    for (int k=gd.kstart+1; k<gd.kend; ++k)
    {
        TF diff = std::abs(gd.z[k] - hhub);
        if (diff < minz)
        {
            minz = diff;
            k_hub = k;
        }
    }

    // Gaussian filter width based on grid spacing
    TF Delta = TF(1.5) * gd.dx; // assume uniform spacing
    TF radius = diam * TF(0.5); // rotor radius

    const int jj = gd.jstride;
    const int kk = gd.kstride;

    // Arrays of grid indices and weights that cover the disk
    indices.clear();
    weights.clear();

    for (int j=gd.jstart; j<gd.jend; ++j)
        for (int i=gd.istart; i<gd.iend; ++i)
        {
            // Distance of cell centre to turbine
            TF dx = gd.x[i] - xpos;
            TF dy = gd.y[j] - ypos;
            TF r2 = dx*dx + dy*dy;
            if (std::sqrt(r2) <= radius)
            {
                // Store flattened index and Gaussian weight
                indices.push_back(i + j*jj + k_hub*kk);
                weights.push_back(std::exp(-6.*r2/(Delta*Delta)));
            }
        }

    // Normalise weights so they sum to one
    TF wsum = std::accumulate(weights.begin(), weights.end(), TF(0));
    if (wsum > TF(0))
        for (auto& w : weights)
            w /= wsum;
}

template<typename TF>
void Turbine<TF>::exec(Stats<TF>&, double time)
{
    // Skip forcing until turbine start time
    if (time < turbstarttime)
        return;

    // Access velocity fields and grid information
    auto& u = fields.mp.at("u")->fld;
    auto& v = fields.mp.at("v")->fld;
    auto& gd = grid.get_grid_data();

    const int jj = gd.jstride;
    const int kk = gd.kstride;

    // Update yaw based on upstream flow measurement
    if (swdynyaw && time >= next_yaw)
    {
        // find cell one diameter upstream of current yaw
        TF xref = xpos - diam * std::cos(yaw);
        TF yref = ypos - diam * std::sin(yaw);

        auto nearest_index = [](TF val, const std::vector<TF>& arr, int s, int e)
        {
            int idx = s;
            TF md = std::abs(arr[s] - val);
            for (int i=s+1; i<e; ++i)
            {
                TF d = std::abs(arr[i] - val);
                if (d < md)
                {
                    md = d;
                    idx = i;
                }
            }
            return idx;
        };

        int iu = nearest_index(xref, gd.x, gd.istart, gd.iend);
        int ju = nearest_index(yref, gd.y, gd.jstart, gd.jend);
        int idx = iu + ju*jj + k_hub*kk;

        TF ur = u[idx];
        TF vr = v[idx];

        TF target = std::atan2(vr, ur);      // target wind direction
        yaw += (target - yaw) * TF(0.2);      // relax toward target

        next_yaw += yawperiod;
    }

    // Calculate disk-averaged incoming velocity
    TF umean = 0;
    for (size_t n=0; n<indices.size(); ++n)
    {
        int idx = indices[n];
        umean += weights[n] * (u[idx]*std::cos(yaw) + v[idx]*std::sin(yaw));
    }

    // Compute thrust force and power from disk-averaged velocity
    TF thrust = 0.5 * ct * umean * std::abs(umean);
    TF powloc = cp * 0.5 * umean * umean * umean * area;
    power = powloc;

    // Apply actuator disk forcing to all covered cells
    for (size_t n=0; n<indices.size(); ++n)
    {
        int idx = indices[n];
        TF f = thrust * weights[n];
        u[idx] -= f * std::cos(yaw);
        v[idx] -= f * std::sin(yaw);
    }
}

template<typename TF>
void Turbine<TF>::prepare_device()
{
#ifdef USECUDA
    // Allocate GPU resources (placeholder)
    dummy_g.allocate(1);
#endif
}

template<typename TF>
void Turbine<TF>::clear_device()
{
#ifdef USECUDA
    // Release GPU resources (placeholder)
    dummy_g.clear();
#endif
}

#ifdef FLOAT_SINGLE
template class Turbine<float>;
#else
template class Turbine<double>;
#endif
