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

#ifndef DRAGDISK_H
#define DRAGDISK_H

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
class DragDisk
{
    public:
        DragDisk(Master&, Grid<TF>&, Fields<TF>&, Input&);
        ~DragDisk();

        void init(Input&);                  ///< Initialize drag disk from input file.
        void create(Input&, Stats<TF>&);    ///< Setup drag disk indices.
        void exec();                        ///< Apply drag forcing.

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
        bool sw_disk;        ///< switch for drag disk
        bool sw_diskstart;   ///< switch for drag disk start 

        // Drag disk settings
        TF diameter;   ///< disk diameter expressed in horizontal grid cells
        TF height;      ///< disk center height in meters
        TF cd;          ///< drag coefficient
        TF xc;          ///< streamwise coordinates of disk center (m)
        TF yc;          ///< horizontal coordinates of disk center (m)
        TF zc;          ///< vertical coordinate of disk center (m)
        TF radius;      ///< disk radius expressed in grid cells
        int k_center;   ///< vertical index of disk center
        int i_center;   ///< horizontal i-index of disk center
        int j_center;   ///< horizontal j-index of disk center
        std::vector<int> disk_indices; ///< flattened indices of points inside disk

        bool has_custom_diameter = false;
        bool has_custom_height = false;
        bool has_custom_xc = false;
        bool has_custom_yc = false;

        #ifdef USECUDA
        #endif
};

#endif
