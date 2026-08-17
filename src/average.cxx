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

#include <cstdio>
#include <fstream>
#include <string>

#include "average.h"
#include "constants.h"
#include "dump.h"
#include "grid.h"
#include "input.h"
#include "master.h"
#include "timeloop.h"

template<typename TF>
Average<TF>::Average(Master& masterin, Grid<TF>& gridin, Dump<TF>& dumpin, Input& inputin) :
    master(masterin), grid(gridin), dump(dumpin)
{
    swaverage = inputin.get_item<bool>("average", "swaverage", "", false);

    if (swaverage)
    {
        starttime = inputin.get_item<double>("average", "starttime", "", 0.);
        sampletime = inputin.get_item<double>("average", "sampletime", "");
        averagelist = inputin.get_list<std::string>(
                "average", "averagelist", "", std::vector<std::string>());

        if (averagelist.empty())
        {
            std::string msg = "Empty Average list";
            throw std::runtime_error(msg);
        }
    }
    else
    {
        inputin.flag_as_used("average", "starttime", "");
        inputin.flag_as_used("average", "sampletime", "");
        inputin.flag_as_used("average", "averagelist", "");
    }
}

template<typename TF>
Average<TF>::~Average()
{
}

template<typename TF>
void Average<TF>::init()
{
    if (!swaverage)
        return;

    istarttime = convert_to_itime(starttime);
    isampletime = convert_to_itime(sampletime);
    ieffective_starttime = istarttime;
    average_time = 0.;
}

template<typename TF>
void Average<TF>::create()
{
    if (!swaverage)
        return;

    if (!averagelist.empty())
    {
        for (auto& it : averagelist)
            master.print_warning("field %s in [average][averagelist] is illegal\n", it.c_str());
    }
}

template<typename TF>
unsigned long Average<TF>::get_time_limit(unsigned long itime)
{
    if (!swaverage || !has_fields())
        return Constants::ulhuge;

    if (itime < ieffective_starttime)
        return ieffective_starttime - itime;

    const unsigned long remainder = (itime - istarttime) % isampletime;

    if (remainder == 0)
        return isampletime;
    else
        return isampletime - remainder;
}

template<typename TF>
bool Average<TF>::has_started(unsigned long itime) const
{
    return itime > ieffective_starttime;
}

template<typename TF>
std::vector<std::string>& Average<TF>::get_averagelist()
{
    return averagelist;
}

template<typename TF>
void Average<TF>::add_field(const std::string& name)
{
    if (!swaverage || average_fields.find(name) != average_fields.end())
        return;

    auto& gd = grid.get_grid_data();
    average_fields[name].fld.resize(gd.ncells, TF(0));
}

template<typename TF>
void Average<TF>::accumulate(const std::string& name, const TF* const data, const TF dt)
{
    if (!swaverage)
        return;

    auto& average_field = average_fields.at(name);
    auto& gd = grid.get_grid_data();

    for (int n=0; n<gd.ncells; ++n)
        average_field.fld[n] += data[n] * dt;
}

template<typename TF>
bool Average<TF>::is_output_time(unsigned long itime) const
{
    return itime > istarttime && ((itime - istarttime) % isampletime == 0);
}

template<typename TF>
bool Average<TF>::has_average_state_file(const std::string& name, int iotime) const
{
    char filename[256];
    std::snprintf(filename, 256, "%s.%07d", name.c_str(), iotime);

    std::ifstream infile(filename);

    return infile.good();
}

template<typename TF>
void Average<TF>::start_fresh_from_restart(unsigned long itime)
{
    reset_fields_host();
    average_time = 0.;

    if (itime > istarttime)
        ieffective_starttime = itime;
    else
        ieffective_starttime = istarttime;
}

template<typename TF>
void Average<TF>::reset_fields_host()
{
    auto& gd = grid.get_grid_data();

    for (auto& it : average_fields)
        std::fill(it.second.fld.begin(), it.second.fld.begin() + gd.ncells, TF(0));
}

template<typename TF>
void Average<TF>::finish_step(const TF dt, unsigned long itime, int iotime)
{
    if (!swaverage || !has_fields())
        return;

    if (itime <= ieffective_starttime)
    {
        reset_fields();
        average_time = 0.;
        return;
    }

    average_time += dt;

    if (!is_output_time(itime))
        return;

    for (auto& it : average_fields)
        write_field(it.first, iotime);

    reset_fields();
    average_time = 0.;
}

#ifndef USECUDA
template<typename TF>
void Average<TF>::save(const int iotime)
{
    if (!swaverage || !has_fields())
        return;

    for (auto& it : average_fields)
        dump.save_data(it.second.fld.data(), it.first + "_avg_acc", iotime);

    save_average_time(iotime);
}
#endif

template<typename TF>
void Average<TF>::load(const int iotime, unsigned long itime)
{
    if (!swaverage || !has_fields())
        return;

    int nfiles = has_average_state_file("average_time", iotime) ? 1 : 0;
    const int ntotal = static_cast<int>(average_fields.size()) + 1;

    for (auto& it : average_fields)
        nfiles += has_average_state_file(it.first + "_avg_acc", iotime) ? 1 : 0;

    if (nfiles == 0)
    {
        start_fresh_from_restart(itime);

        if (iotime > 0)
        {
            master.print_warning(
                    "Average restart state missing at time %07d, starting fresh averages from t=%lu with [average][starttime]=%lu\n",
                    iotime, ieffective_starttime, istarttime);
        }

        return;
    }

    if (nfiles != ntotal)
        throw std::runtime_error("Incomplete average restart state");

    for (auto& it : average_fields)
        dump.load_data(it.second.fld.data(), it.first + "_avg_acc", iotime);

    load_average_time(iotime);
    ieffective_starttime = istarttime;
}

#ifndef USECUDA
template<typename TF>
void Average<TF>::reset_fields()
{
    reset_fields_host();
}

template<typename TF>
void Average<TF>::write_field(const std::string& name, int iotime)
{
    auto& average_field = average_fields.at(name);
    auto& gd = grid.get_grid_data();
    const TF fac = TF(1) / TF(average_time);

    for (int n=0; n<gd.ncells; ++n)
        average_field.fld[n] *= fac;

    dump.save_data(average_field.fld.data(), name + "_avg", iotime);
}
#endif

template<typename TF>
void Average<TF>::save_average_time(const int iotime)
{
    int nerror = 0;

    if (master.get_mpiid() == 0)
    {
        char filename[256];
        std::snprintf(filename, 256, "%s.%07d", "average_time", iotime);

        FILE* pFile = fopen(filename, "wbx");

        if (pFile == NULL)
        {
            ++nerror;
        }
        else
        {
            fwrite(&average_time, sizeof(double), 1, pFile);
            fclose(pFile);
        }
    }

    master.sum(&nerror, 1);

    if (nerror)
        throw std::runtime_error("Error writing average_time");
}

template<typename TF>
void Average<TF>::load_average_time(const int iotime)
{
    int nerror = 0;

    if (master.get_mpiid() == 0)
    {
        char filename[256];
        std::snprintf(filename, 256, "%s.%07d", "average_time", iotime);

        FILE* pFile = fopen(filename, "rb");

        if (pFile == NULL)
        {
            if (iotime == 0)
                average_time = 0.;
            else
                ++nerror;
        }
        else
        {
            if (fread(&average_time, sizeof(double), 1, pFile) != 1u)
                ++nerror;

            fclose(pFile);
        }
    }

    master.sum(&nerror, 1);

    if (nerror)
        throw std::runtime_error("Error reading average_time");

    master.broadcast(&average_time, 1);
}

#ifdef FLOAT_SINGLE
template class Average<float>;
#else
template class Average<double>;
#endif
