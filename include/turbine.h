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
class Netcdf_handle;
template<typename> class Grid;
template<typename> class Fields;
template<typename> class Field3d_operators;
template<typename> class Field3d_io;
template<typename> class Timedep;
template<typename> class Stats;
template<typename> class Thermo;

template<typename TF>
class Turbine
{
    public:
        Turbine(Master&, Grid<TF>&, Fields<TF>&, Input&);
        ~Turbine();

        void init();                        // Initialize turbines from input file.
        void create(Input&, Stats<TF>&);    // Setup disk indicies
        void exec(Stats<TF>& stats);        // Calculate velocities, forcing, power, yaw of actuator disk

        // GPU functions and variables
        #ifdef USECUDA
        void prepare_device();
        void clear_device();
        #endif

    private:
        Master& master;     
        Grid<TF>& grid;     
        Fields<TF>& fields; 
        Field3d_operators<TF> field3d_operators;
        Field3d_io<TF> field3d_io;

        // Internal switches
        bool sw_adm;        // switch to enable turbine modules
        bool sw_windfarm;   // switch to enable multipl turbines in a farm arrangement

        // wind farm parameters
        int nturbine;   // number of turbines if sw_windfarm is enabled

        // turbine parameters
        TF diameter;    // disk diameter [m]
        TF height;      // disk center height in meters
        TF area_disk;   // geometrical disk area [m^2]
        TF cp;          // power coefficient (rotationl forcing)
        TF ct;          // thrust coefficient (axial forcing)
        TF tsr;         // tip-speed ratio
        TF xc;          // streamwise coordinates of disk center [m]
        TF yc;          // horizontal coordinates of disk center [m]
        int k_center;   // vertical index of disk center
        int i_center;   // horizontal i-index of disk center
        int j_center;   // horizontal j-index of disk center
        std::vector<int> turbine_idx;   // Turbine index
        std::vector<int> disk_indices;  // flattened indices of points inside disk
        std::vector<TF>  radial_dist;   // radial distance of points from hub center
        std::vector<TF>  disk_vel;      // space averaged disk velocity [m/s]
        std::vector<TF>  omega;         // radial velocity of rotor [1/s]

        // Turbine statistics/diagnostics
        TF turbine_starttime;   // start time for turbines [s]
        TF turbine_statstart;   // start time for turbine statistics [s]
        TF turbine_statperiod;  // period for turbine statistics [s]

        // Turbine forcing variables
        TF ftx; // thrust force on actuator disk
        TF fty; // tanent force in y, for wake rotation
        TF ftz; // tanent force in z, for wake rotation

        // Strings
        const std::string tend_name = "turbine";    // tendency name
        const std::string layout_file;              // turbine layout file name

        #ifdef USECUDA
        #endif
};

#endif
