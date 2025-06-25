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

#ifndef WINDFARM_H
#define WINDFARM_H

#include <vector>
#include <string>
#include "turbine.h"

class Master;
class Input;
template<typename> class Grid;
template<typename> class Fields;
template<typename> class Stats;

/**
 * Wind farm manager class that configures multiple turbines.
 * Reads the following parameters from (case).ini file
 *
 * [windfarm]
 * nturbrows   ; number of turbine rows
 * nturbcols   ; number of turbine columns
 * spacingx    ; row spacing (x direction) (D)
 * spacingy    ; column spacing (y direction) (D)
 * swstaggered ; enable staggered layout
 * farmlocx    ; farm location x (m)
 * farmlocy    ; farm location y (m)
 */

template<typename TF>
class Windfarm
{
    public:
        Windfarm(Master&, Grid<TF>&, Fields<TF>&, Input&);
        ~Windfarm();

        void create();                ///< Setup wind farm layout.
        void exec(Stats<TF>&, double);///< Apply turbine forcing for all turbines.

        TF get_farm_power() const { return farm_power; }

        void prepare_device();
        void clear_device();

    private:
        Master& master;   ///< Global simulation controller
        Grid<TF>& grid;   ///< LES grid description
        Fields<TF>& fields; ///< Flow field container
        Input& input;     ///< Input file handle

        int nturbrows;    ///< Number of rows in the farm
        int nturbcols;    ///< Number of columns in the farm
        TF spacingx;      ///< Row spacing in diameters
        TF spacingy;      ///< Column spacing in diameters
        bool swstaggered; ///< Whether layout uses staggered rows
        TF farmlocx;      ///< Starting x location of farm (m)
        TF farmlocy;      ///< Starting y location of farm (m)
        std::string layoutfile; ///< Optional manual layout file

        TF diam;          ///< Turbine diameter (m)
        
        TF farm_power;    ///< Aggregate power output (W)

        std::vector<Turbine<TF>> turbines; ///< Container of turbines
};

#endif
