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

#ifndef TURBINE_H
#define TURBINE_H

#include <vector>
#include <string>

class Master;
class Input;
template<typename> class Grid;
template<typename> class Fields;
template<typename> class Stats;

/**
 * Simple actuator disk turbine model.
 * Reads the following parameters from (case).ini file
 *
 * [turbine]
 * diam            ; turbine diameter (m)
 * hhub            ; hub height (m)
 * ct              ; thrust coefficient
 * cp              ; power coefficient
 * tsr             ; tip speed ratio
 * swdynyaw        ; enable dynamic yaw control
 * yawperiod       ; yaw period (s)
 * turbstarttime   ; start time for turbines (s)
 * swturbstats     ; enable turbine statistics
 * turbstatperiod  ; period for turbine statistics (s)
 */

template<typename TF>
class Turbine
{
    public:
        Turbine(Master&, Grid<TF>&, Fields<TF>&, Input&, TF, TF, TF);
        ~Turbine();

        void create();                ///< Setup turbine parameters.
        void exec(Stats<TF>&, double);///< Apply turbine forcing.

        TF get_power() const { return power; }

        // GPU interface
        void prepare_device();
        void clear_device();

    private:
        Master& master;  ///< Reference to global simulation controller
        Grid<TF>& grid;  ///< LES grid description
        Fields<TF>& fields; ///< Flow field container

        TF diam;         ///< Rotor diameter (m)
        TF hhub;         ///< Hub height (m)
        TF ct;           ///< Thrust coefficient
        TF cp;           ///< Power coefficient
        TF tsr;          ///< Tip speed ratio
        bool swdynyaw;   ///< Enable dynamic yaw control
        TF yawperiod;    ///< Yaw update period (s)
        TF turbstarttime;///< Start time for turbine forcing (s)
        bool swturbstats;///< Enable turbine statistics output
        TF turbstatperiod;///< Statistics output period (s)

        // Runtime parameters
        TF xpos;         ///< Turbine x-position (m)
        TF ypos;         ///< Turbine y-position (m)
        TF yaw;          ///< Current yaw angle (rad)
        TF next_yaw;     ///< Next time to update yaw (s)
        TF area;         ///< Disk area (m^2)
        int k_hub;       ///< Grid index of hub height

        std::vector<int> indices; ///< Grid indices covered by disk
        std::vector<TF> weights;  ///< Corresponding Gaussian weights

        TF power;        ///< Instantaneous turbine power (W)

        // GPU containers
        #ifdef USECUDA
        cuda_vector<TF> dummy_g; ///< placeholder for future GPU data
        #endif
};

#endif
