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

#include "windfarm.h"
#include "master.h"
#include "grid.h"
#include "fields.h"
#include "input.h"
#include "stats.h"
#include <fstream>
#include <stdexcept>

template<typename TF>
Windfarm<TF>::Windfarm(Master& masterin, Grid<TF>& gridin, Fields<TF>& fieldsin, Input& inputin) :
    master(masterin), grid(gridin), fields(fieldsin), input(inputin)
{
    // Read wind farm layout parameters
    nturbrows  = inputin.get_item<int>("windfarm", "nturbrows", "", 1);
    nturbcols  = inputin.get_item<int>("windfarm", "nturbcols", "", 1);
    spacingx   = inputin.get_item<TF>("windfarm", "spacingx", "", TF(0));
    spacingy   = inputin.get_item<TF>("windfarm", "spacingy", "", TF(0));
    swstaggered= inputin.get_item<bool>("windfarm", "swstaggered", "", false);
    farmlocx   = inputin.get_item<TF>("windfarm", "farmlocx", "", TF(0));
    farmlocy   = inputin.get_item<TF>("windfarm", "farmlocy", "", TF(0));
    layoutfile = inputin.get_item<std::string>("windfarm", "layoutfile", "", "");
    diam       = inputin.get_item<TF>("turbine", "diam", "");
    farm_power = 0; // initialise aggregate power
}

template<typename TF>
Windfarm<TF>::~Windfarm()
{
}

template<typename TF>
void Windfarm<TF>::create()
{
    // Remove existing turbine list
    turbines.clear();

    auto& gd = grid.get_grid_data();

    // Helper lambda to create a single turbine
    auto add_turb = [&](TF x, TF y, TF hh)
    {
        if (x - 0.5*diam < 0 || x + 0.5*diam > gd.xsize ||
            y - 0.5*diam < 0 || y + 0.5*diam > gd.ysize)
            throw std::runtime_error("Turbine near boundary");

        turbines.emplace_back(master, grid, fields, input, x, y, hh);
        turbines.back().create();
    };

    // Option 1: read turbine positions from external file
    if (!layoutfile.empty())
    {
        std::ifstream f(layoutfile);
        if (!f.is_open())
            throw std::runtime_error("Cannot open layout file " + layoutfile);
        TF x, y, hh;
        while (f >> x >> y >> hh)
            add_turb(x, y, hh);
    }
    else
    {
        // Option 2: generate regular grid of turbines
        for (int r=0; r<nturbrows; ++r)
        {
            TF y = farmlocy + r*spacingy*diam;
            for (int c=0; c<nturbcols; ++c)
            {
                TF offset = (swstaggered && (r%2==1)) ? 0.5*spacingx*diam : TF(0);
                TF x = farmlocx + c*spacingx*diam + offset;
                add_turb(x, y, -1);
            }
        }
    }
}

template<typename TF>
void Windfarm<TF>::exec(Stats<TF>& stats, double time)
{
    farm_power = 0; // reset accumulator
    for (auto& t : turbines)
    {
        t.exec(stats, time);
        farm_power += t.get_power();
    }
}

template<typename TF>
void Windfarm<TF>::prepare_device()
{
    // Allocate GPU data for all turbines
    for (auto& t : turbines)
        t.prepare_device();
}

template<typename TF>
void Windfarm<TF>::clear_device()
{
    // Clear GPU data for all turbines
    for (auto& t : turbines)
        t.clear_device();
}

#ifdef FLOAT_SINGLE
template class Windfarm<float>;
#else
template class Windfarm<double>;
#endif
