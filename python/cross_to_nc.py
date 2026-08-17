#
#  MicroHH
#  Copyright (c) 2011-2024 Chiel van Heerwaarden
#  Copyright (c) 2011-2024 Thijs Heus
#  Copyright (c) 2014-2024 Bart van Stratum
#
#  This file is part of MicroHH
#
#  MicroHH is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  MicroHH is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with MicroHH.  If not, see <http://www.gnu.org/licenses/>.
#

import os
import microhh_tools as mht  # available in microhh/python directory
import argparse
import collections
import glob
import numpy as np
from multiprocessing import Pool, set_start_method
import platform

if platform.system() == 'Darwin':
    set_start_method('fork')

HALF_LEVEL_VARS = {'w', 'lflx', 'sflx'}


def get_source_variable(variable):
    if avg and not variable.endswith('_avg'):
        return '{0}_avg'.format(variable)
    else:
        return variable


def get_base_variable(variable):
    if variable.endswith('_avg'):
        return variable[:-4]
    else:
        return variable


def get_halflevel(variable):
    base_variable = get_base_variable(variable)
    return ''.join([
        '1' if base_variable == 'u' else '0',
        '1' if base_variable == 'v' else '0',
        '1' if base_variable in HALF_LEVEL_VARS else '0'])


def get_vertical_coordinate(variable):
    if get_halflevel(variable)[2] == '1':
        return grid.dim['zh'][:-1]
    else:
        return grid.dim['z']


def get_cross_coordinate(mode, variable):
    halflevel = get_halflevel(variable)

    if mode == 'xy':
        return get_vertical_coordinate(variable)
    elif mode == 'xz':
        return grid.dim['yh'] if halflevel[1] == '1' else grid.dim['y']
    elif mode == 'yz':
        return grid.dim['xh'] if halflevel[0] == '1' else grid.dim['x']
    else:
        raise ValueError('Unsupported mode {}'.format(mode))


def get_average_indexes(mode, variable):
    if indexes is not None:
        return indexes

    values = nl['cross'][mode]
    values = [values] if not isinstance(values, list) else values
    coordinate = get_cross_coordinate(mode, variable)

    return [int(np.abs(coordinate - value).argmin()) for value in values]


def get_first_available_time(variables):
    first_otimes = []

    for variable in variables:
        files = glob.glob('{0}.[0-9]*'.format(get_source_variable(variable)))
        if len(files) == 0:
            raise Exception('Cannot find any average files for {0}'.format(variable))

        otimes = [int(os.path.basename(f).split('.')[-1]) for f in files]
        first_otimes.append(min(otimes))

    return max(first_otimes) * 10**iotimeprec


def _build_time_otime_pairs(starttime, endtime, sampletime, iotimeprec):
    scale = 10**iotimeprec
    times = np.arange(starttime, endtime + sampletime, sampletime)
    return [(time, int(round(time / scale))) for time in times]


def _cross_section_shape(mode, at_surface):
    if at_surface or mode == 'xy':
        return (jtot, itot)
    if mode == 'xz':
        return (ktot, itot)
    return (ktot, jtot)


def _time_block_shape(mode, nindices, at_surface):
    if at_surface:
        return (jtot, itot)
    if mode == 'xy':
        return (nindices, jtot, itot)
    if mode == 'xz':
        return (ktot, nindices, itot)
    return (ktot, jtot, nindices)


def convert_to_nc(variables):
    # Loop over the different variables and crosssections
    for variable in variables:
        source_variable = get_source_variable(variable)
        for mode in modes:
            filename = "{0}.{1}.nc".format(source_variable, mode)
            try:
                otime = time_otime_pairs[0][1]

                if avg:
                    at_surface = False
                    halflevel = get_halflevel(source_variable)
                    indexes_local = get_average_indexes(mode, source_variable)
                else:
                    if os.path.isfile("{0}.xy.000.{1:07d}".format(source_variable, otime)):
                        if mode != 'xy':
                            continue
                        at_surface = True
                    else:
                        at_surface = False

                    halflevel = '000'
                    if not at_surface:
                        if indexes is None:
                            indexes_local,halflevel = mht.get_cross_indices(source_variable, mode)
                        else:
                            indexes_local = indexes

                            files = glob.glob("{0:}.{1}.*.{2:05d}.{3:07d}".format(
                                    source_variable, mode, indexes_local[0], otime))
                            if len(files) == 0:
                                raise Exception('Cannot find any cross-section')
                            halflevel = files[0].split('.')[-3]

                dim = collections.OrderedDict()
                dim['time'] = []
                dim['z'] = range(ktot)
                dim['y'] = range(jtot)
                dim['x'] = range(itot)

                if at_surface:
                    dim.pop('z')
                    n = itot * jtot
                    indexes_local = [-1]
                elif mode == 'xy':
                    dim.update({'z': []})
                    n = itot * jtot
                elif mode == 'xz':
                    dim.update({'y': []})
                    n = itot * ktot
                elif mode == 'yz':
                    dim.update({'x': []})
                    n = ktot * jtot

                section_shape = _cross_section_shape(mode, at_surface)
                block_shape = _time_block_shape(mode, len(indexes_local), at_surface)

                if halflevel[0] == '1':
                    dim['xh'] = dim.pop('x')
                if halflevel[1] == '1':
                    dim['yh'] = dim.pop('y')
                if halflevel[2] == '1':
                    dim['zh'] = dim.pop('z')
                ncfile = mht.Create_ncfile(
                    grid, filename, source_variable, dim, precision, compression)
                
                for key, val in dim.items():
                    if key == 'time':
                        continue
                    elif val == []:
                        ncfile.dimvar[key][:] = grid.dim[key][indexes_local]

                for t, (time, otime) in enumerate(time_otime_pairs):
                    if avg:
                        f_in = "{0:}.{1:07d}".format(source_variable, otime)
                        try:
                            fin = mht.Read_binary(grid, f_in)
                        except Exception as ex:
                            print(ex)
                            continue

                        print("Processing %8s, time=%7i" % (source_variable, otime))

                        field = fin.read(itot * jtot * ktot).reshape((ktot, jtot, itot))
                        time_block = np.empty(block_shape, dtype=grid.dtype)

                        for k, index in enumerate(indexes_local):
                            if mode == 'xy':
                                time_block[k, :, :] = field[index, :, :]
                            elif mode == 'xz':
                                time_block[:, k, :] = field[:, index, :]
                            elif mode == 'yz':
                                time_block[:, :, k] = field[:, :, index]

                        ncfile.dimvar['time'][t] = time
                        ncfile.var[t, ...] = time_block
                        fin.close()
                        continue

                    if at_surface:
                        f_in = "{0}.{1}.{2}.{3:07d}".format(
                            source_variable, mode, halflevel, otime)
                        try:
                            fin = mht.Read_binary(grid, f_in)
                        except Exception as ex:
                            print (ex)
                            continue

                        print(
                            "Processing %8s, time=%7i, index=%4i" %
                            (source_variable, otime, -1))

                        ncfile.dimvar['time'][t] = time
                        ncfile.var[t, :, :] = fin.read(n).reshape(section_shape)
                        fin.close()
                        continue

                    time_block = np.empty(block_shape, dtype=grid.dtype)
                    missing_file = False

                    for k in range(len(indexes_local)):
                        index = indexes_local[k]
                        f_in = "{0:}.{1}.{2}.{3:05d}.{4:07d}".format(
                            source_variable, mode, halflevel, index, otime)
                        try:
                            fin = mht.Read_binary(grid, f_in)
                        except Exception as ex:
                            print (ex)
                            missing_file = True
                            break

                        print(
                            "Processing %8s, time=%7i, index=%4i" %
                            (source_variable, otime, index))

                        data = fin.read(n).reshape(section_shape)
                        if mode == 'xy':
                            time_block[k, :, :] = data
                        elif mode == 'xz':
                            time_block[:, k, :] = data
                        elif mode == 'yz':
                            time_block[:, :, k] = data

                        fin.close()

                    if missing_file:
                        continue

                    ncfile.dimvar['time'][t] = time
                    ncfile.var[t, ...] = time_block
                ncfile.close()

            except Exception as ex:
                print(ex)
                print("Failed to create %s" % filename)



# Parse command line and namelist options
cross_modes = ['xy', 'xz', 'yz']
parser = argparse.ArgumentParser(
    description='Convert MicroHH binary cross-sections to netCDF4 files.')
parser.add_argument(
    '-m',
    '--modes',
    nargs='*',
    help='mode of the cross section',
    choices=cross_modes)
parser.add_argument('-f', '--filename', help='ini file name')
parser.add_argument('-d', '--directory', help='directory')
parser.add_argument('-v', '--vars', nargs='*', help='variable names')
parser.add_argument('-x', '--index', nargs='*', help='indices', type=int)
parser.add_argument('-t0', '--starttime', help='first time step to be parsed')
parser.add_argument('-t1', '--endtime', help='last time step to be parsed')
parser.add_argument(
    '-tstep',
    '--sampletime',
    help='time interval to be parsed')
parser.add_argument(
    '-p',
    '--precision',
    help='precision',
    choices=[
        'single',
         'double'])
parser.add_argument(
    '-n',
    '--nprocs',
    help='Number of processes',
    type=int,
    default=1)
parser.add_argument(
    '-c',
    '--nocompression',
    help='do not compress the netcdf file',
    action='store_true')

parser.add_argument(
    '-o',
    '--order',
    help='order',
    choices=[
        2, 4], type = int)
parser.add_argument(
    '--avg',
    help='target running-average 3D fields',
    action='store_true')

args = parser.parse_args()

if args.directory is not None:
    os.chdir(args.directory)

modes = args.modes
indexes = args.index

nl = mht.Read_namelist(args.filename)
itot = nl['grid']['itot']
jtot = nl['grid']['jtot']
ktot = nl['grid']['ktot']

starttime = float(
    args.starttime) if args.starttime is not None else nl['time']['starttime']
endtime = float(
    args.endtime) if args.endtime is not None else nl['time']['endtime']
avg = args.avg
if avg:
    starttime = float(args.starttime) if args.starttime is not None else nl['average'].get('starttime', 0.)
    sampletime = float(
        args.sampletime) if args.sampletime is not None else nl['average']['sampletime']
else:
    sampletime = float(
        args.sampletime) if args.sampletime is not None else nl['cross']['sampletime']

if args.modes is None:
    modes = list(nl['cross'].keys() & cross_modes)

    # Check if there are paths in the cross-list
    if 'xy' not in modes:
        for v in nl['cross']['crosslist']:
            if 'path' in v:
                modes.append('xy')
                break
else:
    modes = args.modes

if 'iotimeprec' in nl['time']:
    iotimeprec = nl['time']['iotimeprec']
else:
    iotimeprec = 0.

variables = args.vars if args.vars is not None else (
        nl['average']['averagelist'] if avg else nl['cross']['crosslist'])

# In case variables is a single string, convert to list.
variables = [ variables ] if not isinstance(variables, list) else variables

if avg and args.starttime is None:
    starttime = get_first_available_time(variables)

precision = args.precision
nprocs = args.nprocs if args.nprocs is not None else len(variables)
compression = not(args.nocompression)
try:
    order = args.order if args.order is not None else nl['grid']['swspatialorder']
except KeyError:
    order = 2

# End option parsing
grid = mht.Read_grid(itot, jtot, ktot, order = order)
time_otime_pairs = _build_time_otime_pairs(starttime, endtime, sampletime, iotimeprec)

chunks = [variables[i::nprocs] for i in range(nprocs)]

pool = Pool(processes=nprocs)

for _ in pool.imap_unordered(convert_to_nc, chunks):
    pass

pool.close()
pool.join()
