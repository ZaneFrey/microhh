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

 * MicroHH is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.

 * You should have received a copy of the GNU General Public License
 * along with MicroHH.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef DEFINES_H
#define DEFINES_H

#ifdef USECUDA
#include <cuda_runtime.h>
#else
struct dim3
{
    constexpr dim3(const unsigned int xin=1, const unsigned int yin=1, const unsigned int zin=1) :
        x(xin), y(yin), z(zin) {}
    unsigned int x;
    unsigned int y;
    unsigned int z;
};
#endif

#define restrict RESTRICTKEYWORD
enum class Sim_mode { Init, Run, Post };

#endif
