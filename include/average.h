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

#ifndef AVERAGE_H
#define AVERAGE_H

#include <map>
#include <string>
#include <vector>

#include "cuda_buffer.h"

class Master;
class Input;
template<typename> class Grid;
template<typename> class Dump;

template<typename TF>
class Average
{
    public:
        Average(Master&, Grid<TF>&, Dump<TF>&, Input&);
        ~Average();

        void init();
        void create();

        unsigned long get_time_limit(unsigned long);
        bool has_started(unsigned long) const;
        bool get_switch() const { return swaverage; }
        bool has_fields() const { return !average_fields.empty(); }
        std::vector<std::string>& get_averagelist();

        void add_field(const std::string&);
        void accumulate(const std::string&, const TF* const, const TF);
        void finish_step(const TF, unsigned long, int);
        void save(const int);
        void load(const int, unsigned long);

        #ifdef USECUDA
        void prepare_device();
        void clear_device();
        void forward_device();
        void accumulate_g(const std::string&, const TF* const, const TF);
        #endif

    private:
        struct Average_field
        {
            std::vector<TF> fld;

            #ifdef USECUDA
            cuda_vector<TF> fld_g;
            #endif
        };

        Master& master;
        Grid<TF>& grid;
        Dump<TF>& dump;

        std::vector<std::string> averagelist;
        std::map<std::string, Average_field> average_fields;

        bool swaverage;
        double starttime;
        double sampletime;
        unsigned long istarttime;
        unsigned long ieffective_starttime;
        unsigned long isampletime;
        double average_time;

        bool is_output_time(unsigned long) const;
        bool has_average_state_file(const std::string&, int) const;
        void start_fresh_from_restart(unsigned long);
        void reset_fields_host();
        void reset_fields();
        void write_field(const std::string&, int);
        void save_average_time(const int);
        void load_average_time(const int);
};
#endif
