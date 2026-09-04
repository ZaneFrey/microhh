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

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <limits>
#include <set>
#include <sstream>
#include <stdexcept>

#include "constants.h"
#include "fields.h"
#include "grid.h"
#include "input.h"
#include "master.h"
#include "netcdf_interface.h"
#include "timeloop.h"
#include "windfarm.h"
#ifdef USECUDA
#include "windfarm_kernels.cuh"
#endif

namespace
{
    template<typename TF>
    TF wrap_coordinate(TF x, const TF length)
    {
        x -= std::floor(x/length)*length;
        return x == length ? TF(0) : x;
    }

    template<typename TF>
    TF minimum_image(const TF x, const TF length)
    {
        return x - std::round(x/length)*length;
    }

    template<typename TF>
    TF interpolate(const TF* const field, const Windfarm_interpolation<TF>& interpolation)
    {
        TF value = TF(0);
        for (int n=0; n<8; ++n)
            value += interpolation.weight[n]*field[interpolation.index[n]];
        return value;
    }

    void hash_bytes(std::uint64_t& hash, const void* const data, const std::size_t size)
    {
        const unsigned char* bytes = static_cast<const unsigned char*>(data);
        for (std::size_t n=0; n<size; ++n)
        {
            hash ^= bytes[n];
            hash *= UINT64_C(1099511628211);
        }
    }

    template<typename T>
    void hash_value(std::uint64_t& hash, const T value)
    {
        hash_bytes(hash, &value, sizeof(T));
    }

    template<typename TF>
    TF wrap_angle(TF angle)
    {
        angle -= TF(360)*std::floor((angle+TF(180))/TF(360));
        return angle;
    }

    template<typename TF>
    TF angle_difference(const TF target, const TF current)
    {
        return wrap_angle(target-current);
    }
}

template<typename TF>
WindFarm<TF>::WindFarm(
        Master& masterin, Grid<TF>& gridin, Fields<TF>& fieldsin, Input& input, const bool has_dem) :
    master(masterin),
    grid(gridin),
    fields(fieldsin),
    sw_windfarm(Windfarm_type::Disabled),
    diameter(TF(0)),
    hub_height(TF(0)),
    disk_memory_time(TF(600)),
    actuator_resolution(TF(0)),
    start_time(TF(0)),
    ramp_time(TF(0)),
    rotation_sign(TF(1)),
    ct_prime(TF(4)/TF(3)),
    cp_prime(TF(0.8)),
    tip_speed_ratio(TF(8)),
    sample_time(TF(0)),
    yaw_memory_time(TF(60)),
    yaw_update_time(TF(60)),
    yaw_sensor_distance(TF(1)),
    yaw_rate(TF(0.5)),
    yaw_deadband(TF(2)),
    isample_time(0),
    istart_time(0),
    iyaw_update_time(0),
    dynamic_yaw(false),
    configuration_hash(0),
    created(false),
    device_prepared(false),
    suppress_first_output(false),
    sensor_geometry_ready(false),
    mean_warning_printed(false),
    last_update_time(std::numeric_limits<unsigned long>::max()),
    output_record(0)
{
    const std::string sw = input.get_item<std::string>("windfarm", "swwf", "", "0");
    if (sw == "0")
        return;
    else if (sw == "adm")
        sw_windfarm = Windfarm_type::Adm;
    else if (sw == "admr")
        sw_windfarm = Windfarm_type::Admr;
    else
        throw std::runtime_error("Illegal value for [windfarm][swwf]");

    if (has_dem)
        throw std::runtime_error("The wind-farm model is incompatible with DEM immersed boundaries");

    diameter = input.get_item<TF>("windfarm", "diameter", "m");
    hub_height = input.get_item<TF>("windfarm", "hubheight", "m");
    disk_memory_time = input.get_item<TF>("windfarm", "diskmemtime", "s", TF(600));
    actuator_resolution = input.get_item<TF>("windfarm", "adresolution", "m");
    turbine_file = input.get_item<std::string>("windfarm", "wfturblocs", "");
    turbine_format = input.get_item<std::string>("windfarm", "wfturblocformat", "");
    start_time = input.get_item<TF>("windfarm", "wfstarttime", "s", TF(0));
    ramp_time = input.get_item<TF>("windfarm", "wframptime", "s", TF(0));
    yaw_file = input.get_item<std::string>("windfarm", "wfyawfile", "", "");
    sample_time = input.get_item<TF>("windfarm", "sampletime", "s", TF(0));
    ct_prime = input.get_item<TF>("windfarm", "ctprime", "", TF(4)/TF(3));
    dynamic_yaw = input.get_item<bool>("windfarm", "swdynamicyaw", "", false);

    if (dynamic_yaw)
    {
        yaw_memory_time = input.get_item<TF>("windfarm", "yawmemtime", "s", TF(60));
        yaw_update_time = input.get_item<TF>("windfarm", "yawupdatetime", "s", TF(60));
        yaw_sensor_distance = input.get_item<TF>("windfarm", "yawsensordist", "", TF(1));
        yaw_rate = input.get_item<TF>("windfarm", "yawrate", "degree s-1", TF(0.5));
        yaw_deadband = input.get_item<TF>("windfarm", "yawdeadband", "degree", TF(2));
    }

    if (sw_windfarm == Windfarm_type::Admr)
    {
        cp_prime = input.get_item<TF>("windfarm", "cpprime", "", TF(0.8));
        tip_speed_ratio = input.get_item<TF>("windfarm", "tsr", "", TF(8));
        rotation_sign = input.get_item<TF>("windfarm", "rotationsign", "", TF(1));
    }
}

template<typename TF>
WindFarm<TF>::~WindFarm() = default;

template<typename TF>
void WindFarm<TF>::validate_configuration() const
{
    const auto finite = [](const TF value) { return std::isfinite(value); };
    if (!finite(diameter) || diameter <= TF(0) || !finite(hub_height) || hub_height <= TF(0))
        throw std::runtime_error("Wind-turbine diameter and hub height must be finite and positive");
    if (!finite(actuator_resolution) || actuator_resolution <= TF(0))
        throw std::runtime_error("[windfarm][adresolution] must be finite and positive");
    if (!finite(disk_memory_time) || disk_memory_time < TF(0) ||
        !finite(start_time) || start_time < TF(0) ||
        !finite(ramp_time) || ramp_time < TF(0) ||
        !finite(sample_time) || sample_time < TF(0))
        throw std::runtime_error("Wind-farm times must be finite and non-negative");
    if (!finite(ct_prime) || ct_prime < TF(0))
        throw std::runtime_error("[windfarm][ctprime] must be finite and non-negative");
    if (turbine_format != "xy" && turbine_format != "latlong")
        throw std::runtime_error("[windfarm][wfturblocformat] must be either xy or latlong");
    if (turbine_format == "latlong")
    {
        const auto& gd = grid.get_grid_data();
        if (!std::isfinite(gd.lat) || !std::isfinite(gd.lon) ||
            gd.lat < TF(-90) || gd.lat > TF(90) || gd.lon < TF(-180) || gd.lon > TF(180))
            throw std::runtime_error("The latlong turbine format requires valid [grid] lat and lon values");
    }
    if (sw_windfarm == Windfarm_type::Admr &&
        (!finite(cp_prime) || cp_prime < TF(0) || !finite(tip_speed_ratio) || tip_speed_ratio <= TF(0) ||
         (rotation_sign != TF(-1) && rotation_sign != TF(1))))
        throw std::runtime_error("Invalid ADMR power, tip-speed-ratio, or rotation-sign setting");
    if (dynamic_yaw &&
        (!finite(yaw_memory_time) || yaw_memory_time < TF(0) ||
         !finite(yaw_update_time) || yaw_update_time < TF(0) ||
         !finite(yaw_sensor_distance) || yaw_sensor_distance <= TF(0) ||
         !finite(yaw_rate) || yaw_rate <= TF(0) ||
         !finite(yaw_deadband) || yaw_deadband < TF(0) || yaw_deadband >= TF(180)))
        throw std::runtime_error("Invalid dynamic-yaw memory, update, sensor-distance, rate, or deadband setting");
}

template<typename TF>
void WindFarm<TF>::read_turbines()
{
    std::vector<int> ids;
    std::vector<TF> coordinates;
    int read_ok = 1;

    if (master.get_mpiid() == 0)
    {
        try
        {
            std::ifstream file(turbine_file);
            if (!file)
                throw std::runtime_error("Cannot open wind-turbine location file \"" + turbine_file + "\"");

            std::set<int> unique_ids;
            std::string line;
            int line_number = 0;
            while (std::getline(file, line))
            {
                ++line_number;
                line = line.substr(0, line.find('#'));
                if (line.find_first_not_of(" \t\r\n") == std::string::npos)
                    continue;

                int id;
                double first;
                double second;
                std::string extra;
                std::istringstream stream(line);
                if (!(stream >> id >> first >> second) || stream >> extra ||
                    !std::isfinite(first) || !std::isfinite(second))
                    throw std::runtime_error("Invalid turbine-location row " + std::to_string(line_number));
                if (!unique_ids.insert(id).second)
                    throw std::runtime_error("Duplicate wind-turbine ID " + std::to_string(id));

                TF x;
                TF y;
                if (turbine_format == "xy")
                {
                    x = TF(first);
                    y = TF(second);
                }
                else
                {
                    const auto& gd = grid.get_grid_data();
                    const double dlon = second-double(gd.lon);
                    if (std::abs(dlon) > 180.)
                        throw std::runtime_error("A turbine longitude requires unsupported date-line wrapping");
                    const double earth_radius = 6.371e6;
                    x = gd.xsize/TF(2)+TF(earth_radius*std::cos(double(gd.lat)*M_PI/180.)*dlon*M_PI/180.);
                    y = gd.ysize/TF(2)+TF(earth_radius*(first-double(gd.lat))*M_PI/180.);
                }
                ids.push_back(id);
                coordinates.push_back(x);
                coordinates.push_back(y);
            }
            if (ids.empty())
                throw std::runtime_error("The wind-turbine location file is empty");
        }
        catch (const std::exception& exception)
        {
            master.print_warning("%s\n", exception.what());
            read_ok = 0;
        }
    }

    master.broadcast(&read_ok, 1);
    if (!read_ok)
        throw std::runtime_error("Failed to read the wind-turbine location file");

    int count = int(ids.size());
    master.broadcast(&count, 1);
    ids.resize(count);
    coordinates.resize(2*count);
    master.broadcast(ids.data(), count);
    master.broadcast(coordinates.data(), 2*count);

    turbines.resize(count);
    for (int n=0; n<count; ++n)
    {
        Turbine turbine{};
        turbine.id = ids[n];
        turbine.x = coordinates[2*n];
        turbine.y = coordinates[2*n+1];
        turbine.nx = TF(1);
        turbine.last_sensor_sample_time = std::numeric_limits<unsigned long>::max();
        turbine.last_yaw_update_time = std::numeric_limits<unsigned long>::max();
        turbine.last_rotor_sample_time = std::numeric_limits<unsigned long>::max();
        turbines[n] = turbine;
    }
}

template<typename TF>
void WindFarm<TF>::read_yaw()
{
    if (yaw_file.empty())
        return;

    std::vector<int> ids;
    std::vector<TF> yaw;
    int read_ok = 1;
    if (master.get_mpiid() == 0)
    {
        try
        {
            std::ifstream file(yaw_file);
            if (!file)
                throw std::runtime_error("Cannot open wind-turbine yaw file \"" + yaw_file + "\"");
            std::set<int> unique_ids;
            std::string line;
            int line_number = 0;
            while (std::getline(file, line))
            {
                ++line_number;
                line = line.substr(0, line.find('#'));
                if (line.find_first_not_of(" \t\r\n") == std::string::npos)
                    continue;
                int id;
                double angle;
                std::string extra;
                std::istringstream stream(line);
                if (!(stream >> id >> angle) || stream >> extra || !std::isfinite(angle))
                    throw std::runtime_error("Invalid turbine-yaw row " + std::to_string(line_number));
                if (!unique_ids.insert(id).second)
                    throw std::runtime_error("Duplicate wind-turbine yaw ID " + std::to_string(id));
                ids.push_back(id);
                yaw.push_back(TF(angle));
            }
            if (ids.size() != turbines.size())
                throw std::runtime_error("The yaw file must contain exactly one row per turbine");
        }
        catch (const std::exception& exception)
        {
            master.print_warning("%s\n", exception.what());
            read_ok = 0;
        }
    }
    master.broadcast(&read_ok, 1);
    if (!read_ok)
        throw std::runtime_error("Failed to read the wind-turbine yaw file");

    int count = int(ids.size());
    master.broadcast(&count, 1);
    ids.resize(count);
    yaw.resize(count);
    master.broadcast(ids.data(), count);
    master.broadcast(yaw.data(), count);

    std::set<int> matched;
    for (int n=0; n<count; ++n)
    {
        auto turbine = std::find_if(turbines.begin(), turbines.end(),
                [id=ids[n]](const Turbine& item) { return item.id == id; });
        if (turbine == turbines.end() || !matched.insert(ids[n]).second)
            throw std::runtime_error("The yaw file contains an unknown or duplicate turbine ID");
        turbine->yaw_offset = wrap_angle(yaw[n]);
    }
    if (matched.size() != turbines.size())
        throw std::runtime_error("The yaw file does not list every turbine ID");

    for (auto& turbine : turbines)
        turbine.rotor_yaw = turbine.yaw_offset;
}

template<typename TF>
void WindFarm<TF>::build_topology()
{
    const auto& gd = grid.get_grid_data();
    const TF radius = diameter/TF(2);
    if (diameter >= gd.xsize || diameter >= gd.ysize)
        throw std::runtime_error("A wind-turbine rotor overlaps its own periodic image");
    if (hub_height-radius <= TF(0) || hub_height+radius >= gd.zsize)
        throw std::runtime_error("A wind-turbine rotor intersects a vertical domain boundary");

    for (auto& turbine : turbines)
    {
        if (!(turbine.x > TF(0) && turbine.x < gd.xsize && turbine.y > TF(0) && turbine.y < gd.ysize))
            throw std::runtime_error("Wind-turbine hubs must lie strictly inside the horizontal domain");
        turbine.x = wrap_coordinate(turbine.x, gd.xsize);
        turbine.y = wrap_coordinate(turbine.y, gd.ysize);
    }

    for (std::size_t a=0; a<turbines.size(); ++a)
        for (std::size_t b=a+1; b<turbines.size(); ++b)
        {
            const TF dx = minimum_image(turbines[a].x-turbines[b].x, gd.xsize);
            const TF dy = minimum_image(turbines[a].y-turbines[b].y, gd.ysize);
            if (std::sqrt(dx*dx+dy*dy) < diameter)
                throw std::runtime_error("Wind-turbine rotor disks overlap under periodic minimum-image distance");
        }

    elements.clear();
    for (std::size_t t=0; t<turbines.size(); ++t)
    {
        Turbine& turbine = turbines[t];
        turbine.element_start = int(elements.size());
        const int radial_count = std::max(1, int(std::round(radius/actuator_resolution)));
        const TF dr = radius/TF(radial_count);
        for (int j=0; j<radial_count; ++j)
        {
            const TF r = (TF(j)+TF(0.5))*dr;
            const int azimuth_count = std::max(18, int(std::round(TF(2*M_PI)*r/actuator_resolution)));
            const TF dtheta = TF(2*M_PI)/TF(azimuth_count);
            const TF area = r*dr*dtheta;
            for (int l=0; l<azimuth_count; ++l)
            {
                const TF theta = (TF(l)+TF(0.5))*dtheta;
                Windfarm_element<TF> element{};
                element.turbine = int(t);
                element.ring = j;
                element.point = l;
                element.r = r;
                element.theta = theta;
                element.area = area;
                elements.push_back(element);
            }
        }
        turbine.element_count = int(elements.size())-turbine.element_start;
        TF actuator_area = TF(0);
        for (int n=0; n<turbine.element_count; ++n)
            actuator_area += elements[turbine.element_start+n].area;
        const TF rotor_area = TF(M_PI)*radius*radius;
        if (std::abs(actuator_area-rotor_area) >
            TF(100)*std::numeric_limits<TF>::epsilon()*rotor_area)
            throw std::runtime_error("The discrete actuator area does not equal the physical rotor area");
    }
    element_forces.resize(elements.size());
}

template<typename TF>
void WindFarm<TF>::update_rotor_geometry()
{
    const auto& gd = grid.get_grid_data();
    const TF radius = diameter/TF(2);
    for (auto& turbine : turbines)
    {
        turbine.rotor_yaw = wrap_angle(turbine.rotor_yaw);
        const TF yaw_radians = turbine.rotor_yaw*TF(M_PI/180.);
        turbine.nx = std::cos(yaw_radians);
        turbine.ny = std::sin(yaw_radians);
        const TF epsilon_n = std::sqrt(std::pow(turbine.nx*gd.dx, 2)+std::pow(turbine.ny*gd.dy, 2));
        const TF epsilon_1 = std::sqrt(std::pow(turbine.ny*gd.dx, 2)+std::pow(turbine.nx*gd.dy, 2));
        const TF gaussian_x = TF(3)*(std::abs(turbine.nx)*epsilon_n+std::abs(turbine.ny)*epsilon_1);
        const TF gaussian_y = TF(3)*(std::abs(turbine.ny)*epsilon_n+std::abs(turbine.nx)*epsilon_1);
        if (TF(2)*(radius*std::abs(turbine.ny)+gaussian_x) >= gd.xsize ||
            TF(2)*(radius*std::abs(turbine.nx)+gaussian_y) >= gd.ysize)
            throw std::runtime_error("A turbine rotor and Gaussian support overlap their own periodic image");
    }

    for (auto& element : elements)
    {
        const Turbine& turbine = turbines[element.turbine];
        const TF e1x = -turbine.ny;
        const TF e1y = turbine.nx;
        const TF er_x = std::cos(element.theta)*e1x;
        const TF er_y = std::cos(element.theta)*e1y;
        const TF er_z = std::sin(element.theta);
        element.x = wrap_coordinate(turbine.x+element.r*er_x, gd.xsize);
        element.y = wrap_coordinate(turbine.y+element.r*er_y, gd.ysize);
        element.z = hub_height+element.r*er_z;
        element.etheta_x = turbine.ny*er_z;
        element.etheta_y = -turbine.nx*er_z;
        element.etheta_z = turbine.nx*er_y-turbine.ny*er_x;
    }
}

template<typename TF>
Windfarm_interpolation<TF> WindFarm<TF>::make_interpolation(
        const TF x, const TF y, const TF z, const int x_face, const int y_face, const int z_face) const
{
    const auto& gd = grid.get_grid_data();
    const TF sx = x/gd.dx-(x_face ? TF(0) : TF(0.5));
    const TF sy = y/gd.dy-(y_face ? TF(0) : TF(0.5));
    const int i0_global = int(std::floor(sx));
    const int j0_global = int(std::floor(sy));
    const TF wx = sx-TF(i0_global);
    const TF wy = sy-TF(j0_global);
    const int global_i_start = master.get_MPI_data().mpicoordx*gd.imax;
    const int global_j_start = master.get_MPI_data().mpicoordy*gd.jmax;
    const int i0 = i0_global-global_i_start+gd.igc;
    const int j0 = j0_global-global_j_start+gd.jgc;

    const std::vector<TF>& vertical = z_face ? gd.zh : gd.z;
    int k0 = gd.kstart-1;
    while (k0+1 < gd.kend+1 && vertical[k0+1] <= z)
        ++k0;
    if (k0 < 0 || k0+1 >= gd.kcells || vertical[k0+1] == vertical[k0])
        throw std::runtime_error("Cannot bracket an actuator point on the vertical grid");
    const TF wz = (z-vertical[k0])/(vertical[k0+1]-vertical[k0]);

    Windfarm_interpolation<TF> result;
    int n = 0;
    for (int kk=0; kk<2; ++kk)
        for (int jj=0; jj<2; ++jj)
            for (int ii=0; ii<2; ++ii)
            {
                result.index[n] = (i0+ii)+(j0+jj)*gd.icells+(k0+kk)*gd.ijcells;
                result.weight[n] = (ii ? wx : TF(1)-wx)*(jj ? wy : TF(1)-wy)*(kk ? wz : TF(1)-wz);
                ++n;
            }
    return result;
}

template<typename TF>
TF WindFarm<TF>::interpolate_reference_density(const TF height) const
{
    const auto& gd = grid.get_grid_data();
    int k0 = gd.kstart-1;
    while (k0+1 < gd.kend+1 && gd.z[k0+1] <= height)
        ++k0;
    const TF weight = (height-gd.z[k0])/(gd.z[k0+1]-gd.z[k0]);
    return (TF(1)-weight)*fields.rhoref[k0]+weight*fields.rhoref[k0+1];
}

template<typename TF>
TF WindFarm<TF>::local_dz(const TF height) const
{
    const auto& gd = grid.get_grid_data();
    int k0 = gd.kstart-1;
    while (k0+1 < gd.kend+1 && gd.z[k0+1] <= height)
        ++k0;
    const TF weight = (height-gd.z[k0])/(gd.z[k0+1]-gd.z[k0]);
    return (TF(1)-weight)*gd.dz[k0]+weight*gd.dz[k0+1];
}

template<typename TF>
void WindFarm<TF>::build_sampling()
{
    const auto& gd = grid.get_grid_data();
    samples.clear();
    for (const auto& element : elements)
    {
        int global_i = std::min(gd.itot-1, int(std::floor(element.x/gd.dx)));
        int global_j = std::min(gd.jtot-1, int(std::floor(element.y/gd.dy)));
        const int owner_x = global_i/gd.imax;
        const int owner_y = global_j/gd.jmax;
        Windfarm_sample<TF> sample;
        sample.owner = owner_x == master.get_MPI_data().mpicoordx && owner_y == master.get_MPI_data().mpicoordy;
        sample.turbine = element.turbine;
        sample.area = element.area;
        sample.density = interpolate_reference_density(element.z);
        sample.nx = turbines[element.turbine].nx;
        sample.ny = turbines[element.turbine].ny;
        if (sample.owner)
        {
            sample.u = make_interpolation(element.x, element.y, element.z, 1, 0, 0);
            sample.v = make_interpolation(element.x, element.y, element.z, 0, 1, 0);
            sample.w = make_interpolation(element.x, element.y, element.z, 0, 0, 1);
        }
        samples.push_back(sample);
    }
    sums.resize(3*turbines.size());
}

template<typename TF>
void WindFarm<TF>::build_sensor_sampling()
{
    const auto& gd = grid.get_grid_data();
    sensor_samples.clear();
    sensor_samples.reserve(elements.size());
    for (const auto& element : elements)
    {
        const Turbine& turbine = turbines[element.turbine];
        const TF yaw_radians = turbine.yaw_command*TF(M_PI/180.);
        const TF nx = std::cos(yaw_radians);
        const TF ny = std::sin(yaw_radians);
        const TF center_x = turbine.x-yaw_sensor_distance*diameter*nx;
        const TF center_y = turbine.y-yaw_sensor_distance*diameter*ny;
        const TF er_x = -std::cos(element.theta)*ny;
        const TF er_y = std::cos(element.theta)*nx;
        const TF er_z = std::sin(element.theta);
        const TF x = center_x+element.r*er_x;
        const TF y = center_y+element.r*er_y;
        const TF z = hub_height+element.r*er_z;
        if (!(x > TF(0) && x < gd.xsize && y > TF(0) && y < gd.ysize))
            throw std::runtime_error(
                    "Dynamic-yaw sensor disk for turbine " + std::to_string(turbine.id) +
                    " crosses a horizontal domain boundary");

        const int global_i = std::min(gd.itot-1, int(std::floor(x/gd.dx)));
        const int global_j = std::min(gd.jtot-1, int(std::floor(y/gd.dy)));
        const int owner_x = global_i/gd.imax;
        const int owner_y = global_j/gd.jmax;
        Windfarm_sensor_sample<TF> sample{};
        sample.owner = owner_x == master.get_MPI_data().mpicoordx &&
                       owner_y == master.get_MPI_data().mpicoordy;
        sample.turbine = element.turbine;
        sample.area = element.area;
        if (sample.owner)
        {
            sample.u = make_interpolation(x, y, z, 1, 0, 0);
            sample.v = make_interpolation(x, y, z, 0, 1, 0);
        }
        sensor_samples.push_back(sample);
    }
    sensor_sums.resize(3*turbines.size());
    sensor_geometry_ready = true;
}

template<typename TF>
void WindFarm<TF>::build_scatter()
{
    const auto& gd = grid.get_grid_data();
    const int global_i_start = master.get_MPI_data().mpicoordx*gd.imax;
    const int global_j_start = master.get_MPI_data().mpicoordy*gd.jmax;
    std::vector<TF> normalization(3*elements.size(), TF(0));
    for (auto& component : scatter)
        component.clear();

    for (std::size_t p=0; p<elements.size(); ++p)
    {
        const auto& element = elements[p];
        const Turbine& turbine = turbines[element.turbine];
        const TF dz = local_dz(element.z);
        const TF epsilon_n = std::sqrt(std::pow(turbine.nx*gd.dx, 2)+std::pow(turbine.ny*gd.dy, 2));
        const TF epsilon_1 = std::sqrt(std::pow(turbine.ny*gd.dx, 2)+std::pow(turbine.nx*gd.dy, 2));
        const TF epsilon_2 = dz;
        const TF hx = TF(3)*(std::abs(turbine.nx)*epsilon_n+std::abs(turbine.ny)*epsilon_1);
        const TF hy = TF(3)*(std::abs(turbine.ny)*epsilon_n+std::abs(turbine.nx)*epsilon_1);
        const TF hz = TF(3)*epsilon_2;
        if (TF(2)*hx >= gd.xsize || TF(2)*hy >= gd.ysize)
            throw std::runtime_error("A turbine Gaussian stencil overlaps its own periodic image");
        if (element.z-hz <= TF(0) || element.z+hz >= gd.zsize)
            throw std::runtime_error("A turbine Gaussian stencil intersects a vertical domain boundary");

        for (int component=0; component<3; ++component)
        {
            const TF x_offset = component == 0 ? TF(0) : TF(0.5);
            const TF y_offset = component == 1 ? TF(0) : TF(0.5);
            const int i_first = int(std::ceil((element.x-hx)/gd.dx-x_offset));
            const int i_last = int(std::floor((element.x+hx)/gd.dx-x_offset));
            const int j_first = int(std::ceil((element.y-hy)/gd.dy-y_offset));
            const int j_last = int(std::floor((element.y+hy)/gd.dy-y_offset));
            const int k_first = component == 2 ? gd.kstart+1 : gd.kstart;
            const int k_last = gd.kend;
            const std::vector<TF>& vertical = component == 2 ? gd.zh : gd.z;

            for (int k=k_first; k<k_last; ++k)
            {
                if (vertical[k] < element.z-hz || vertical[k] > element.z+hz)
                    continue;
                for (int j_unwrapped=j_first; j_unwrapped<=j_last; ++j_unwrapped)
                    for (int i_unwrapped=i_first; i_unwrapped<=i_last; ++i_unwrapped)
                    {
                        const int global_i = (i_unwrapped%gd.itot+gd.itot)%gd.itot;
                        const int global_j = (j_unwrapped%gd.jtot+gd.jtot)%gd.jtot;
                        if (global_i < global_i_start || global_i >= global_i_start+gd.imax ||
                            global_j < global_j_start || global_j >= global_j_start+gd.jmax)
                            continue;
                        const TF x = (TF(i_unwrapped)+x_offset)*gd.dx;
                        const TF y = (TF(j_unwrapped)+y_offset)*gd.dy;
                        const TF dx = minimum_image(x-element.x, gd.xsize);
                        const TF dy = minimum_image(y-element.y, gd.ysize);
                        const TF dzp = vertical[k]-element.z;
                        const TF dn = dx*turbine.nx+dy*turbine.ny;
                        const TF d1 = -dx*turbine.ny+dy*turbine.nx;
                        if (std::abs(dn) > TF(3)*epsilon_n || std::abs(d1) > TF(3)*epsilon_1 ||
                            std::abs(dzp) > TF(3)*epsilon_2)
                            continue;
                        const TF kernel = std::exp(TF(-0.5)*(dn*dn/(epsilon_n*epsilon_n)+
                                                               d1*d1/(epsilon_1*epsilon_1)+
                                                               dzp*dzp/(epsilon_2*epsilon_2)));
                        const TF volume = gd.dx*gd.dy*(component == 2 ? gd.dzh[k] : gd.dz[k]);
                        normalization[3*p+component] += kernel*volume;
                        const int i = global_i-global_i_start+gd.igc;
                        const int j = global_j-global_j_start+gd.jgc;
                        const int index = i+j*gd.icells+k*gd.ijcells;
                        const TF density = component == 2 ? fields.rhorefh[k] : fields.rhoref[k];
                        if (!(density > TF(0)) || !std::isfinite(density))
                            throw std::runtime_error("Reference density must be finite and positive in a turbine stencil");
                        scatter[component].push_back({int(p), index, kernel/density});
                    }
            }
        }
    }

    master.sum(normalization.data(), int(normalization.size()));
    for (int component=0; component<3; ++component)
        for (auto& entry : scatter[component])
        {
            const TF norm = normalization[3*entry.element+component];
            if (!(norm > TF(0)) || !std::isfinite(norm))
                throw std::runtime_error("A wind-turbine Gaussian stencil has zero normalization");
            entry.coefficient /= norm;
        }
}

template<typename TF>
void WindFarm<TF>::calculate_configuration_hash()
{
    configuration_hash = UINT64_C(1469598103934665603);
    const int model = int(sw_windfarm);
    hash_value(configuration_hash, model);
    hash_value(configuration_hash, dynamic_yaw);
    for (const TF value : {diameter, hub_height, disk_memory_time, actuator_resolution,
                            start_time, ramp_time, rotation_sign, ct_prime, cp_prime, tip_speed_ratio,
                            yaw_memory_time, yaw_update_time, yaw_sensor_distance, yaw_rate, yaw_deadband})
    {
        const float compact = float(value);
        hash_value(configuration_hash, compact);
    }
    for (const auto& turbine : turbines)
    {
        hash_value(configuration_hash, turbine.id);
        for (const TF value : {turbine.x, turbine.y, turbine.yaw_offset})
        {
            const float compact = float(value);
            hash_value(configuration_hash, compact);
        }
    }
}

template<typename TF>
void WindFarm<TF>::init()
{
    if (!get_switch())
        return;
    const int ncells = grid.get_grid_data().ncells;
    for (auto& component : source)
        component.assign(ncells, TF(0));
}

template<typename TF>
void WindFarm<TF>::create(Timeloop<TF>& timeloop, const bool is_run)
{
    if (!get_switch() || !is_run)
        return;
    validate_configuration();
    read_turbines();
    read_yaw();
    build_topology();
    calculate_configuration_hash();

    isample_time = sample_time > TF(0) ? convert_to_itime(sample_time) : 0;
    istart_time = convert_to_itime(start_time);
    iyaw_update_time = yaw_update_time > TF(0) ? convert_to_itime(yaw_update_time) : 0;
    if (dynamic_yaw && yaw_update_time > TF(0) && iyaw_update_time == 0)
        throw std::runtime_error("[windfarm][yawupdatetime] is smaller than the time precision");
    load_restart(timeloop.get_iotime(), timeloop.get_itime());
    update_rotor_geometry();
    build_sampling();
    build_scatter();
    if (dynamic_yaw && turbines.front().last_sensor_sample_time != std::numeric_limits<unsigned long>::max())
        build_sensor_sampling();
    suppress_first_output = suppress_first_output && is_sample_time(timeloop.get_itime());
    create_output(timeloop.get_iotime());
    created = true;
    master.print_message("Created wind farm with %d turbines and %d actuator points\n",
                         int(turbines.size()), int(elements.size()));
}

template<typename TF>
void WindFarm<TF>::update_rotor_host()
{
    std::fill(sums.begin(), sums.end(), TF(0));
    const TF* const u = fields.mp.at("u")->fld.data();
    const TF* const v = fields.mp.at("v")->fld.data();
    const TF* const w = fields.mp.at("w")->fld.data();
    const auto& gd = grid.get_grid_data();
    for (std::size_t p=0; p<samples.size(); ++p)
    {
        const auto& sample = samples[p];
        if (!sample.owner)
            continue;
        const TF velocity = (interpolate(u, sample.u)+gd.utrans)*sample.nx+
                            (interpolate(v, sample.v)+gd.vtrans)*sample.ny+
                            TF(0)*interpolate(w, sample.w);
        sums[3*sample.turbine] += velocity*sample.area;
        sums[3*sample.turbine+1] += sample.density*sample.area;
        sums[3*sample.turbine+2] += sample.area;
    }
    master.sum(sums.data(), int(sums.size()));
}

template<typename TF>
void WindFarm<TF>::update_sensor_host()
{
    std::fill(sensor_sums.begin(), sensor_sums.end(), TF(0));
    const TF* const u = fields.mp.at("u")->fld.data();
    const TF* const v = fields.mp.at("v")->fld.data();
    const auto& gd = grid.get_grid_data();
    for (const auto& sample : sensor_samples)
    {
        if (!sample.owner)
            continue;
        sensor_sums[3*sample.turbine] += (interpolate(u, sample.u)+gd.utrans)*sample.area;
        sensor_sums[3*sample.turbine+1] += (interpolate(v, sample.v)+gd.vtrans)*sample.area;
        sensor_sums[3*sample.turbine+2] += sample.area;
    }
    master.sum(sensor_sums.data(), int(sensor_sums.size()));
}

template<typename TF>
void WindFarm<TF>::calculate_hub_mean_host(TF& u_mean, TF& v_mean)
{
    const auto& gd = grid.get_grid_data();
    int k0 = gd.kstart-1;
    while (k0+1 < gd.kend+1 && gd.z[k0+1] <= hub_height)
        ++k0;
    if (k0 < 0 || k0+1 >= gd.kcells || gd.z[k0+1] == gd.z[k0])
        throw std::runtime_error("Cannot bracket the turbine hub height on the vertical grid");
    const TF weight = (hub_height-gd.z[k0])/(gd.z[k0+1]-gd.z[k0]);
    const TF* const u = fields.mp.at("u")->fld.data();
    const TF* const v = fields.mp.at("v")->fld.data();
    hub_sums.assign(2, TF(0));
    for (int j=gd.jstart; j<gd.jend; ++j)
        for (int i=gd.istart; i<gd.iend; ++i)
        {
            const int index0 = i+j*gd.icells+k0*gd.ijcells;
            const int index1 = index0+gd.ijcells;
            hub_sums[0] += (TF(1)-weight)*u[index0]+weight*u[index1]+gd.utrans;
            hub_sums[1] += (TF(1)-weight)*v[index0]+weight*v[index1]+gd.vtrans;
        }
    master.sum(hub_sums.data(), int(hub_sums.size()));
    const TF horizontal_count = TF(gd.itot)*TF(gd.jtot);
    u_mean = hub_sums[0]/horizontal_count;
    v_mean = hub_sums[1]/horizontal_count;
}

template<typename TF>
void WindFarm<TF>::calculate_forces(const TF time, const unsigned long itime)
{
    const TF disk_area = TF(M_PI)*diameter*diameter/TF(4);
    TF gamma = TF(0);
    if (time >= start_time)
        gamma = ramp_time > TF(0) ? TF(1)-std::exp(-(time-start_time)/ramp_time) : TF(1);
    gamma = std::max(TF(0), std::min(TF(1), gamma));
    for (std::size_t t=0; t<turbines.size(); ++t)
    {
        Turbine& turbine = turbines[t];
        const TF actual_area = sums[3*t+2];
        if (!(actual_area > TF(0)))
            throw std::runtime_error("A wind turbine has no owned actuator sampling area");
        turbine.raw_velocity = sums[3*t]/actual_area;
        turbine.density = sums[3*t+1]/actual_area;
        if (!turbine.filter_valid)
        {
            turbine.filtered_velocity = turbine.raw_velocity;
            turbine.filter_valid = true;
        }
        else
        {
            const TF elapsed = TF(itime-turbine.last_rotor_sample_time)/TF(ifactor);
            const TF alpha = disk_memory_time > TF(0) ?
                    TF(1)-std::exp(-elapsed/disk_memory_time) : TF(1);
            turbine.filtered_velocity += alpha*(turbine.raw_velocity-turbine.filtered_velocity);
        }
        turbine.last_rotor_sample_time = itime;

        const TF load_velocity = turbine.filtered_velocity > TF(1.e-6) ? turbine.filtered_velocity : TF(0);
        turbine.thrust = gamma*TF(0.5)*turbine.density*disk_area*ct_prime*load_velocity*load_velocity;
        if (sw_windfarm == Windfarm_type::Admr && load_velocity > TF(0))
        {
            turbine.angular_velocity = tip_speed_ratio*load_velocity/(diameter/TF(2));
            turbine.power = TF(0.5)*turbine.density*disk_area*cp_prime*
                            load_velocity*load_velocity*load_velocity;
            turbine.torque = gamma*turbine.power/turbine.angular_velocity;
        }
        else
            turbine.angular_velocity = turbine.power = turbine.torque = TF(0);

        TF discrete_torque = TF(0);
        for (int n=0; n<turbine.element_count; ++n)
        {
            const auto& element = elements[turbine.element_start+n];
            discrete_torque += element.r*(turbine.torque/actual_area)*(element.area/element.r);
        }
        const TF torque_scale = discrete_torque > TF(0) ? turbine.torque/discrete_torque : TF(0);
        for (int n=0; n<turbine.element_count; ++n)
        {
            const int p = turbine.element_start+n;
            const auto& element = elements[p];
            const TF axial = turbine.thrust*element.area/actual_area;
            const TF tangential = torque_scale*(turbine.torque/actual_area)*(element.area/element.r);
            element_forces[p].x = -axial*turbine.nx-rotation_sign*tangential*element.etheta_x;
            element_forces[p].y = -axial*turbine.ny-rotation_sign*tangential*element.etheta_y;
            element_forces[p].z = -rotation_sign*tangential*element.etheta_z;
        }
    }

    for (auto& component : source)
        std::fill(component.begin(), component.end(), TF(0));
    for (int component=0; component<3; ++component)
        for (const auto& entry : scatter[component])
        {
            const Windfarm_force<TF>& force = element_forces[entry.element];
            const TF value = component == 0 ? force.x : component == 1 ? force.y : force.z;
            source[component][entry.index] += value*entry.coefficient;
        }
}

template<typename TF>
void WindFarm<TF>::initialize_dynamic_yaw(const unsigned long itime)
{
    TF u_mean;
    TF v_mean;
    #ifdef USECUDA
    calculate_hub_mean_g(u_mean, v_mean);
    #else
    calculate_hub_mean_host(u_mean, v_mean);
    #endif
    const TF speed = std::sqrt(u_mean*u_mean+v_mean*v_mean);
    TF initial_command = TF(0);
    if (speed > TF(1.e-6))
        initial_command = wrap_angle(std::atan2(v_mean, u_mean)*TF(180./M_PI));
    else if (!mean_warning_printed)
    {
        master.print_warning(
                "Horizontal mean wind at turbine hub height is too small to initialize dynamic yaw; "
                "retaining a zero-degree command\n");
        mean_warning_printed = true;
    }
    for (auto& turbine : turbines)
        turbine.yaw_command = initial_command;
    build_sensor_sampling();
    #ifdef USECUDA
    refresh_sensor_device();
    update_sensor_g();
    #else
    update_sensor_host();
    #endif
    filter_sensor(itime);
    update_yaw_controller(itime);
}

template<typename TF>
void WindFarm<TF>::filter_sensor(const unsigned long itime)
{
    for (std::size_t t=0; t<turbines.size(); ++t)
    {
        Turbine& turbine = turbines[t];
        const TF actual_area = sensor_sums[3*t+2];
        if (!(actual_area > TF(0)) || !std::isfinite(actual_area))
            throw std::runtime_error("A dynamic-yaw sensor disk has no owned sampling area");
        turbine.sensor_u_raw = sensor_sums[3*t]/actual_area;
        turbine.sensor_v_raw = sensor_sums[3*t+1]/actual_area;
        if (!std::isfinite(turbine.sensor_u_raw) || !std::isfinite(turbine.sensor_v_raw))
            throw std::runtime_error("A dynamic-yaw sensor disk produced a non-finite velocity");
        if (!turbine.sensor_filter_valid)
        {
            turbine.sensor_u_filtered = turbine.sensor_u_raw;
            turbine.sensor_v_filtered = turbine.sensor_v_raw;
            turbine.sensor_filter_valid = true;
        }
        else
        {
            const TF elapsed = TF(itime-turbine.last_sensor_sample_time)/TF(ifactor);
            const TF alpha = yaw_memory_time > TF(0) ?
                    TF(1)-std::exp(-elapsed/yaw_memory_time) : TF(1);
            turbine.sensor_u_filtered += alpha*(turbine.sensor_u_raw-turbine.sensor_u_filtered);
            turbine.sensor_v_filtered += alpha*(turbine.sensor_v_raw-turbine.sensor_v_filtered);
        }
        turbine.last_sensor_sample_time = itime;
    }
}

template<typename TF>
void WindFarm<TF>::update_yaw_controller(const unsigned long itime)
{
    std::vector<TF> commands(turbines.size());
    std::vector<TF> rotor_yaw(turbines.size());
    std::vector<signed char> initialized(turbines.size());
    for (std::size_t t=0; t<turbines.size(); ++t)
    {
        commands[t] = turbines[t].yaw_command;
        rotor_yaw[t] = turbines[t].rotor_yaw;
        initialized[t] = turbines[t].yaw_initialized;
    }

    if (master.get_mpiid() == 0)
        for (std::size_t t=0; t<turbines.size(); ++t)
        {
            const Turbine& turbine = turbines[t];
            const TF u = turbine.sensor_u_filtered;
            const TF v = turbine.sensor_v_filtered;
            const TF speed = std::sqrt(u*u+v*v);
            if (!turbine.sensor_filter_valid || speed <= TF(1.e-6))
                continue;
            const TF command_radians = turbine.yaw_command*TF(M_PI/180.);
            if (u*std::cos(command_radians)+v*std::sin(command_radians) <= TF(0))
                continue;

            const TF candidate = wrap_angle(std::atan2(v, u)*TF(180./M_PI));
            const TF desired = wrap_angle(candidate+turbine.yaw_offset);
            const TF error = angle_difference(desired, turbine.rotor_yaw);
            if (turbine.yaw_initialized && std::abs(error) <= yaw_deadband)
                continue;

            commands[t] = candidate;
            if (!turbine.yaw_initialized)
            {
                rotor_yaw[t] = desired;
                initialized[t] = 1;
            }
            else
            {
                const TF elapsed = TF(itime-turbine.last_yaw_update_time)/TF(ifactor);
                const TF limit = yaw_rate*elapsed;
                rotor_yaw[t] = wrap_angle(turbine.rotor_yaw+std::max(-limit, std::min(limit, error)));
            }
        }

    master.broadcast(commands.data(), int(commands.size()));
    master.broadcast(rotor_yaw.data(), int(rotor_yaw.size()));
    master.broadcast(initialized.data(), int(initialized.size()));
    bool sensor_changed = false;
    bool rotor_changed = false;
    for (std::size_t t=0; t<turbines.size(); ++t)
    {
        sensor_changed = sensor_changed || commands[t] != turbines[t].yaw_command;
        rotor_changed = rotor_changed || rotor_yaw[t] != turbines[t].rotor_yaw;
        turbines[t].yaw_command = wrap_angle(commands[t]);
        turbines[t].rotor_yaw = wrap_angle(rotor_yaw[t]);
        turbines[t].yaw_initialized = initialized[t];
        turbines[t].last_yaw_update_time = itime;
    }
    rebuild_yaw_geometry(sensor_changed, rotor_changed);
}

template<typename TF>
void WindFarm<TF>::rebuild_yaw_geometry(const bool sensor_changed, const bool rotor_changed)
{
    if (sensor_changed)
    {
        build_sensor_sampling();
        #ifdef USECUDA
        refresh_sensor_device();
        #endif
    }
    if (rotor_changed)
    {
        update_rotor_geometry();
        build_sampling();
        build_scatter();
        #ifdef USECUDA
        refresh_rotor_device();
        #endif
    }
}

template<typename TF>
bool WindFarm<TF>::is_sample_time(const unsigned long itime) const
{
    return itime >= istart_time && (isample_time == 0 || (itime-istart_time)%isample_time == 0);
}

template<typename TF>
bool WindFarm<TF>::is_yaw_update_time(const unsigned long itime) const
{
    return dynamic_yaw && itime >= istart_time &&
           (iyaw_update_time == 0 || (itime-istart_time)%iyaw_update_time == 0);
}

template<typename TF>
void WindFarm<TF>::update(Timeloop<TF>& timeloop)
{
    if (!created || timeloop.in_substep() || timeloop.get_itime() == last_update_time)
        return;
    const unsigned long itime = timeloop.get_itime();
    if (dynamic_yaw && itime >= istart_time)
    {
        if (!sensor_geometry_ready)
            initialize_dynamic_yaw(itime);
        else
        {
            #ifdef USECUDA
            update_sensor_g();
            #else
            update_sensor_host();
            #endif
            filter_sensor(itime);
            if (is_yaw_update_time(itime))
                update_yaw_controller(itime);
        }
    }
    #ifdef USECUDA
    update_rotor_g();
    #else
    update_rotor_host();
    #endif
    calculate_forces(TF(timeloop.get_time()), itime);
    #ifdef USECUDA
    cuda_copy(element_forces, element_forces_g);
    for (int component=0; component<3; ++component)
        launch_windfarm_scatter_g(source_g[component].data(), scatter_g[component].data(),
                                int(scatter_g[component].size()), element_forces_g.data(),
                                grid.get_grid_data().ncells, component);
    #endif
    last_update_time = itime;
    if (is_sample_time(itime))
    {
        if (suppress_first_output)
            suppress_first_output = false;
        else
            write_output(timeloop);
    }
}

#ifndef USECUDA
template<typename TF>
void WindFarm<TF>::exec()
{
    if (!created)
        return;
    const auto& gd = grid.get_grid_data();
    TF* const ut = fields.mt.at("u")->fld.data();
    TF* const vt = fields.mt.at("v")->fld.data();
    TF* const wt = fields.mt.at("w")->fld.data();
    for (int k=gd.kstart; k<gd.kend; ++k)
        for (int j=gd.jstart; j<gd.jend; ++j)
            for (int i=gd.istart; i<gd.iend; ++i)
            {
                const int index = i+j*gd.icells+k*gd.ijcells;
                ut[index] += source[0][index];
                vt[index] += source[1][index];
                if (k > gd.kstart)
                    wt[index] += source[2][index];
            }
}
#endif

template<typename TF>
unsigned long WindFarm<TF>::get_time_limit(const unsigned long itime) const
{
    if (!created)
        return Constants::ulhuge;
    if (itime < istart_time)
        return istart_time-itime;
    unsigned long limit = Constants::ulhuge;
    if (isample_time > 0)
    {
        const unsigned long remainder = (itime-istart_time)%isample_time;
        limit = remainder == 0 ? isample_time : isample_time-remainder;
    }
    if (dynamic_yaw && iyaw_update_time > 0)
    {
        const unsigned long remainder = (itime-istart_time)%iyaw_update_time;
        const unsigned long yaw_limit = remainder == 0 ?
                iyaw_update_time : iyaw_update_time-remainder;
        limit = std::min(limit, yaw_limit);
    }
    return limit;
}

template<typename TF>
void WindFarm<TF>::create_output(const int iotime)
{
    char filename[64];
    std::snprintf(filename, sizeof(filename), "windfarm.%07d.nc", iotime);
    output_file = std::make_unique<Netcdf_file>(master, filename, Netcdf_mode::Create);
    output_file->add_dimension("time");
    output_file->add_dimension("turbine", int(turbines.size()));
    iter_var = std::make_unique<Netcdf_variable<int>>(output_file->template add_variable<int>("iter", {"time"}));
    time_var = std::make_unique<Netcdf_variable<TF>>(output_file->template add_variable<TF>("time", {"time"}));
    id_var = std::make_unique<Netcdf_variable<int>>(output_file->template add_variable<int>("id", {"turbine"}));
    yaw_var = std::make_unique<Netcdf_variable<TF>>(output_file->template add_variable<TF>("yaw", {"time", "turbine"}));
    yaw_command_var = std::make_unique<Netcdf_variable<TF>>(
            output_file->template add_variable<TF>("yaw_command", {"time", "turbine"}));
    yaw_offset_var = std::make_unique<Netcdf_variable<TF>>(
            output_file->template add_variable<TF>("yaw_offset", {"turbine"}));
    sensor_u_var = std::make_unique<Netcdf_variable<TF>>(
            output_file->template add_variable<TF>("sensor_u", {"time", "turbine"}));
    sensor_v_var = std::make_unique<Netcdf_variable<TF>>(
            output_file->template add_variable<TF>("sensor_v", {"time", "turbine"}));
    raw_velocity_var = std::make_unique<Netcdf_variable<TF>>(output_file->template add_variable<TF>("raw_velocity", {"time", "turbine"}));
    filtered_velocity_var = std::make_unique<Netcdf_variable<TF>>(output_file->template add_variable<TF>("filtered_velocity", {"time", "turbine"}));
    thrust_var = std::make_unique<Netcdf_variable<TF>>(output_file->template add_variable<TF>("thrust", {"time", "turbine"}));
    torque_var = std::make_unique<Netcdf_variable<TF>>(output_file->template add_variable<TF>("torque", {"time", "turbine"}));
    power_var = std::make_unique<Netcdf_variable<TF>>(output_file->template add_variable<TF>("power", {"time", "turbine"}));

    time_var->add_attribute("units", "s");
    time_var->add_attribute("long_name", "Physical time");
    yaw_var->add_attribute("units", "degree");
    yaw_var->add_attribute("long_name", "Physical turbine rotor yaw angle");
    yaw_command_var->add_attribute("units", "degree");
    yaw_command_var->add_attribute("long_name", "Dynamic-yaw commanded wind direction");
    yaw_offset_var->add_attribute("units", "degree");
    yaw_offset_var->add_attribute("long_name", "Prescribed turbine yaw offset");
    sensor_u_var->add_attribute("units", "m s-1");
    sensor_u_var->add_attribute("long_name", "Filtered upstream-sensor u velocity");
    sensor_v_var->add_attribute("units", "m s-1");
    sensor_v_var->add_attribute("long_name", "Filtered upstream-sensor v velocity");
    raw_velocity_var->add_attribute("units", "m s-1");
    raw_velocity_var->add_attribute("long_name", "Raw disk-normal velocity");
    filtered_velocity_var->add_attribute("units", "m s-1");
    filtered_velocity_var->add_attribute("long_name", "Filtered disk-normal velocity");
    thrust_var->add_attribute("units", "N");
    thrust_var->add_attribute("long_name", "Applied rotor thrust");
    torque_var->add_attribute("units", "N m");
    torque_var->add_attribute("long_name", "Applied rotor torque");
    power_var->add_attribute("units", "W");
    power_var->add_attribute("long_name", "Applied rotor power");
    std::vector<int> ids(turbines.size());
    std::vector<TF> yaw_offset(turbines.size());
    for (std::size_t n=0; n<turbines.size(); ++n)
    {
        ids[n] = turbines[n].id;
        yaw_offset[n] = turbines[n].yaw_offset;
    }
    id_var->insert(ids, {0}, {int(ids.size())});
    yaw_offset_var->insert(yaw_offset, {0}, {int(yaw_offset.size())});
}

template<typename TF>
void WindFarm<TF>::write_output(const Timeloop<TF>& timeloop)
{
    std::vector<TF> yaw(turbines.size());
    std::vector<TF> yaw_command(turbines.size());
    std::vector<TF> sensor_u(turbines.size());
    std::vector<TF> sensor_v(turbines.size());
    std::vector<TF> raw(turbines.size());
    std::vector<TF> filtered(turbines.size());
    std::vector<TF> thrust(turbines.size());
    std::vector<TF> torque(turbines.size());
    std::vector<TF> power(turbines.size());
    for (std::size_t n=0; n<turbines.size(); ++n)
    {
        yaw[n] = turbines[n].rotor_yaw;
        yaw_command[n] = turbines[n].yaw_command;
        sensor_u[n] = turbines[n].sensor_u_filtered;
        sensor_v[n] = turbines[n].sensor_v_filtered;
        raw[n] = turbines[n].raw_velocity;
        filtered[n] = turbines[n].filtered_velocity;
        thrust[n] = turbines[n].thrust;
        torque[n] = turbines[n].torque;
        power[n] = turbines[n].power;
    }
    const int record = int(output_record++);
    const std::vector<int> start = {record, 0};
    const std::vector<int> count = {1, int(turbines.size())};
    iter_var->insert(timeloop.get_iteration(), {record});
    time_var->insert(TF(timeloop.get_time()), {record});
    yaw_var->insert(yaw, start, count);
    yaw_command_var->insert(yaw_command, start, count);
    sensor_u_var->insert(sensor_u, start, count);
    sensor_v_var->insert(sensor_v, start, count);
    raw_velocity_var->insert(raw, start, count);
    filtered_velocity_var->insert(filtered, start, count);
    thrust_var->insert(thrust, start, count);
    torque_var->insert(torque, start, count);
    power_var->insert(power, start, count);
    output_file->sync();
}

template<typename TF>
void WindFarm<TF>::save(const int iotime)
{
    if (!created)
        return;
    if (master.get_mpiid() != 0)
        return;
    char filename[64];
    std::snprintf(filename, sizeof(filename), "windfarm_restart.%07d", iotime);
    std::ofstream file(filename, std::ios::binary | std::ios::trunc);
    if (!file)
        throw std::runtime_error("Cannot create wind-farm restart file");
    const char magic[8] = {'M','H','H','W','F','R','S','T'};
    const std::uint32_t version = 2;
    const std::uint32_t count = std::uint32_t(turbines.size());
    file.write(magic, sizeof(magic));
    file.write(reinterpret_cast<const char*>(&version), sizeof(version));
    file.write(reinterpret_cast<const char*>(&count), sizeof(count));
    file.write(reinterpret_cast<const char*>(&configuration_hash), sizeof(configuration_hash));
    for (const auto& turbine : turbines)
    {
        const std::int32_t id = turbine.id;
        const float yaw_offset = float(turbine.yaw_offset);
        const float yaw_command = float(turbine.yaw_command);
        const float rotor_yaw = float(turbine.rotor_yaw);
        const float sensor_u = float(turbine.sensor_u_filtered);
        const float sensor_v = float(turbine.sensor_v_filtered);
        const std::uint8_t sensor_valid = turbine.sensor_filter_valid;
        const std::uint8_t yaw_initialized = turbine.yaw_initialized;
        const float velocity = float(turbine.filtered_velocity);
        const std::uint8_t valid = turbine.filter_valid;
        const std::uint64_t last_sensor = turbine.last_sensor_sample_time;
        const std::uint64_t last_yaw = turbine.last_yaw_update_time;
        const std::uint64_t last_rotor = turbine.last_rotor_sample_time;
        file.write(reinterpret_cast<const char*>(&id), sizeof(id));
        file.write(reinterpret_cast<const char*>(&yaw_offset), sizeof(yaw_offset));
        file.write(reinterpret_cast<const char*>(&yaw_command), sizeof(yaw_command));
        file.write(reinterpret_cast<const char*>(&rotor_yaw), sizeof(rotor_yaw));
        file.write(reinterpret_cast<const char*>(&sensor_u), sizeof(sensor_u));
        file.write(reinterpret_cast<const char*>(&sensor_v), sizeof(sensor_v));
        file.write(reinterpret_cast<const char*>(&sensor_valid), sizeof(sensor_valid));
        file.write(reinterpret_cast<const char*>(&yaw_initialized), sizeof(yaw_initialized));
        file.write(reinterpret_cast<const char*>(&velocity), sizeof(velocity));
        file.write(reinterpret_cast<const char*>(&valid), sizeof(valid));
        file.write(reinterpret_cast<const char*>(&last_sensor), sizeof(last_sensor));
        file.write(reinterpret_cast<const char*>(&last_yaw), sizeof(last_yaw));
        file.write(reinterpret_cast<const char*>(&last_rotor), sizeof(last_rotor));
    }
    if (!file)
        throw std::runtime_error("Failed while writing the wind-farm restart file");
}

template<typename TF>
void WindFarm<TF>::load_restart(const int iotime, const unsigned long itime)
{
    if (iotime == 0)
        return;
    std::vector<float> state(6*turbines.size());
    std::vector<signed char> flags(3*turbines.size());
    std::vector<unsigned long> times(3*turbines.size());
    int read_ok = 1;
    if (master.get_mpiid() == 0)
    {
        char filename[64];
        std::snprintf(filename, sizeof(filename), "windfarm_restart.%07d", iotime);
        std::ifstream file(filename, std::ios::binary);
        char magic[8];
        std::uint32_t version;
        std::uint32_t count;
        std::uint64_t hash;
        if (!file || !file.read(magic, sizeof(magic)) || !file.read(reinterpret_cast<char*>(&version), sizeof(version)) ||
            !file.read(reinterpret_cast<char*>(&count), sizeof(count)) ||
            !file.read(reinterpret_cast<char*>(&hash), sizeof(hash)) ||
            std::memcmp(magic, "MHHWFRST", 8) != 0 || version != 2 || count != turbines.size() ||
            hash != configuration_hash)
            read_ok = 0;
        for (std::size_t n=0; read_ok && n<turbines.size(); ++n)
        {
            std::int32_t id;
            float yaw_offset;
            std::uint8_t sensor_valid;
            std::uint8_t yaw_initialized;
            std::uint8_t is_valid;
            std::uint64_t last_sensor;
            std::uint64_t last_yaw;
            std::uint64_t last_rotor;
            if (!file.read(reinterpret_cast<char*>(&id), sizeof(id)) ||
                !file.read(reinterpret_cast<char*>(&yaw_offset), sizeof(yaw_offset)) ||
                !file.read(reinterpret_cast<char*>(&state[6*n+1]), sizeof(float)) ||
                !file.read(reinterpret_cast<char*>(&state[6*n+2]), sizeof(float)) ||
                !file.read(reinterpret_cast<char*>(&state[6*n+3]), sizeof(float)) ||
                !file.read(reinterpret_cast<char*>(&state[6*n+4]), sizeof(float)) ||
                !file.read(reinterpret_cast<char*>(&sensor_valid), sizeof(sensor_valid)) ||
                !file.read(reinterpret_cast<char*>(&yaw_initialized), sizeof(yaw_initialized)) ||
                !file.read(reinterpret_cast<char*>(&state[6*n+5]), sizeof(float)) ||
                !file.read(reinterpret_cast<char*>(&is_valid), sizeof(is_valid)) ||
                !file.read(reinterpret_cast<char*>(&last_sensor), sizeof(last_sensor)) ||
                !file.read(reinterpret_cast<char*>(&last_yaw), sizeof(last_yaw)) ||
                !file.read(reinterpret_cast<char*>(&last_rotor), sizeof(last_rotor)) ||
                id != turbines[n].id || yaw_offset != float(turbines[n].yaw_offset) ||
                !std::isfinite(yaw_offset) || !std::isfinite(state[6*n+1]) ||
                !std::isfinite(state[6*n+2]) || !std::isfinite(state[6*n+3]) ||
                !std::isfinite(state[6*n+4]) || !std::isfinite(state[6*n+5]) ||
                state[6*n+1] < -180.f || state[6*n+1] >= 180.f ||
                state[6*n+2] < -180.f || state[6*n+2] >= 180.f ||
                sensor_valid > 1 || yaw_initialized > 1 || is_valid > 1 ||
                (sensor_valid != (last_sensor != std::numeric_limits<std::uint64_t>::max())) ||
                (is_valid != (last_rotor != std::numeric_limits<std::uint64_t>::max())) ||
                (sensor_valid && (last_sensor < istart_time || last_sensor > itime)) ||
                (last_yaw != std::numeric_limits<std::uint64_t>::max() &&
                 (last_yaw < istart_time || last_yaw > itime)) ||
                (yaw_initialized && last_yaw == std::numeric_limits<std::uint64_t>::max()) ||
                (is_valid && last_rotor > itime) || yaw_initialized > sensor_valid ||
                (!dynamic_yaw && (sensor_valid || yaw_initialized ||
                                  last_yaw != std::numeric_limits<std::uint64_t>::max() ||
                                  state[6*n+1] != 0.f || state[6*n+2] != yaw_offset)))
                read_ok = 0;
            else
            {
                state[6*n] = yaw_offset;
                flags[3*n] = static_cast<signed char>(sensor_valid);
                flags[3*n+1] = static_cast<signed char>(yaw_initialized);
                flags[3*n+2] = static_cast<signed char>(is_valid);
                times[3*n] = static_cast<unsigned long>(last_sensor);
                times[3*n+1] = static_cast<unsigned long>(last_yaw);
                times[3*n+2] = static_cast<unsigned long>(last_rotor);
            }
        }
    }
    master.broadcast(&read_ok, 1);
    if (!read_ok)
        throw std::runtime_error("Missing, corrupt, or incompatible wind-farm restart state");
    master.broadcast(state.data(), int(state.size()));
    master.broadcast(flags.data(), int(flags.size()));
    master.broadcast(times.data(), int(times.size()));
    for (std::size_t n=0; n<turbines.size(); ++n)
    {
        turbines[n].yaw_command = TF(state[6*n+1]);
        turbines[n].rotor_yaw = TF(state[6*n+2]);
        turbines[n].sensor_u_filtered = TF(state[6*n+3]);
        turbines[n].sensor_v_filtered = TF(state[6*n+4]);
        turbines[n].sensor_filter_valid = flags[3*n];
        turbines[n].yaw_initialized = flags[3*n+1];
        turbines[n].filtered_velocity = TF(state[6*n+5]);
        turbines[n].filter_valid = flags[3*n+2];
        turbines[n].last_sensor_sample_time = times[3*n];
        turbines[n].last_yaw_update_time = times[3*n+1];
        turbines[n].last_rotor_sample_time = times[3*n+2];
    }
    suppress_first_output = true;
}

#ifdef FLOAT_SINGLE
template class WindFarm<float>;
#else
template class WindFarm<double>;
#endif
