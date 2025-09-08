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
#include "cuda_buffer.h"

class Master;
class Input;
template<typename> class Grid;
template<typename> class Fields;
template<typename> class Stats;

/**
 * Simple circular drag disk that applies a linear damping force
 * to the streamwise velocity inside a specified region.
 *
 * Parameters read from (case).ini:
 *
 * [dragdisk]
 * diameter ; disk diameter in number of grid cells
 * height   ; height of disk centre (m)
 * cd       ; drag coefficient (s^-1, optional, default 0.2)
 */

template<typename TF>
class DragDisk
{
    public:
        DragDisk(Master&, Grid<TF>&, Fields<TF>&, Input&);
        ~DragDisk();

        void create();                ///< Setup drag disk indices.
        void exec(Stats<TF>&, double);///< Apply drag forcing.

        // GPU interface (placeholders)
        void prepare_device();
        void clear_device();

    private:
        Master& master;  ///< Reference to master controller
        Grid<TF>& grid;  ///< Computational grid
        Fields<TF>& fields; ///< Flow fields

        int diameter; ///< disk diameter in grid cells
        TF height;   ///< disk centre height
        TF cd;       ///< drag coefficient
        bool enabled;///< switch for drag disk

        int k_center;             ///< vertical index of disk
        std::vector<int> indices; ///< flattened grid indices inside disk

        #ifdef USECUDA
        cuda_vector<int> indices_g; ///< indices of drag cells on device
        #endif
};

#endif
