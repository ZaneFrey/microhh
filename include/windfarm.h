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

#ifndef WINDFARM_H
#define WINDFARM_H

#include <array>
#include <cstdint>
#include <memory>
#include <string>
#include <vector>

#include "cuda_buffer.h"

class Input;
class Master;
class Netcdf_file;
template<typename> class Fields;
template<typename> class Grid;
template<typename> class Netcdf_variable;
template<typename> class Timeloop;

enum class Windfarm_type { Disabled, Adm, Admr };

template<typename TF>
struct Windfarm_interpolation
{
    int index[8];
    TF weight[8];
};

template<typename TF>
struct Windfarm_sample
{
    int owner;
    int turbine;
    TF area;
    TF density;
    TF nx;
    TF ny;
    Windfarm_interpolation<TF> u;
    Windfarm_interpolation<TF> v;
    Windfarm_interpolation<TF> w;
};

template<typename TF>
struct Windfarm_element
{
    int turbine;
    int ring;
    int point;
    TF x;
    TF y;
    TF z;
    TF r;
    TF theta;
    TF area;
    TF etheta_x;
    TF etheta_y;
    TF etheta_z;
};

template<typename TF>
struct Windfarm_scatter
{
    int element;
    int index;
    TF coefficient;
};

template<typename TF>
struct Windfarm_force
{
    TF x;
    TF y;
    TF z;
};

template<typename TF>
class WindFarm
{
    public:
        WindFarm(Master&, Grid<TF>&, Fields<TF>&, Input&, bool);
        ~WindFarm();

        void init();
        void create(Timeloop<TF>&, bool);
        void update(Timeloop<TF>&);
        void exec();
        void save(int);

        unsigned long get_time_limit(unsigned long) const;
        bool get_switch() const { return sw_windfarm != Windfarm_type::Disabled; }

        #ifdef USECUDA
        void prepare_device();
        void clear_device();
        #endif

    private:
        struct Turbine
        {
            int id;
            int element_start;
            int element_count;
            TF x;
            TF y;
            TF yaw;
            TF nx;
            TF ny;
            TF raw_velocity;
            TF filtered_velocity;
            TF density;
            TF thrust;
            TF torque;
            TF power;
            TF angular_velocity;
            bool filter_valid;
        };

        Master& master;
        Grid<TF>& grid;
        Fields<TF>& fields;

        Windfarm_type sw_windfarm;
        std::string turbine_file;
        std::string turbine_format;
        std::string yaw_file;

        TF diameter;
        TF hub_height;
        TF disk_memory_time;
        TF actuator_resolution;
        TF start_time;
        TF ramp_time;
        TF rotation_sign;
        TF ct_prime;
        TF cp_prime;
        TF tip_speed_ratio;
        TF sample_time;
        unsigned long isample_time;
        unsigned long istart_time;

        std::vector<Turbine> turbines;
        std::vector<Windfarm_element<TF>> elements;
        std::vector<Windfarm_sample<TF>> samples;
        std::array<std::vector<Windfarm_scatter<TF>>, 3> scatter;
        std::vector<Windfarm_force<TF>> element_forces;
        std::array<std::vector<TF>, 3> source;
        std::vector<TF> sums;

        std::uint64_t configuration_hash;
        bool created;
        bool device_prepared;
        bool suppress_first_output;
        unsigned long last_update_time;
        unsigned long output_record;

        std::unique_ptr<Netcdf_file> output_file;
        std::unique_ptr<Netcdf_variable<int>> iter_var;
        std::unique_ptr<Netcdf_variable<TF>> time_var;
        std::unique_ptr<Netcdf_variable<int>> id_var;
        std::unique_ptr<Netcdf_variable<TF>> yaw_var;
        std::unique_ptr<Netcdf_variable<TF>> raw_velocity_var;
        std::unique_ptr<Netcdf_variable<TF>> filtered_velocity_var;
        std::unique_ptr<Netcdf_variable<TF>> thrust_var;
        std::unique_ptr<Netcdf_variable<TF>> torque_var;
        std::unique_ptr<Netcdf_variable<TF>> power_var;

        #ifdef USECUDA
        cuda_vector<Windfarm_sample<TF>> samples_g;
        std::array<cuda_vector<Windfarm_scatter<TF>>, 3> scatter_g;
        cuda_vector<Windfarm_force<TF>> element_forces_g;
        std::array<cuda_vector<TF>, 3> source_g;
        cuda_vector<TF> sums_g;

        void update_g();
        void exec_g();
        #endif

        void read_turbines();
        void read_yaw();
        void validate_configuration() const;
        void build_geometry();
        void build_sampling();
        void build_scatter();
        void calculate_configuration_hash();
        void load_restart(int);
        void create_output(int);
        void write_output(const Timeloop<TF>&);
        void calculate_forces(TF, TF);
        void update_host();

        Windfarm_interpolation<TF> make_interpolation(TF, TF, TF, int, int, int) const;
        TF interpolate_reference_density(TF) const;
        TF local_dz(TF) const;
        bool is_sample_time(unsigned long) const;
};

#endif
