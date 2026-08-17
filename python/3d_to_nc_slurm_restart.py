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

import argparse
import collections
import glob
import os
import platform
import struct
import sys
from dataclasses import dataclass
from multiprocessing import Pool, set_start_method
from pathlib import Path
from typing import Dict, List, Optional

import microhh_tools as mht  # available in microhh/python directory
import numpy as np

if platform.system() == 'Darwin':
    set_start_method('fork')


HALF_LEVEL_VARS = {'w', 'lflx', 'sflx'}


def get_source_variable(variable: str) -> str:
    if avg and not variable.endswith('_avg'):
        return '{0}_avg'.format(variable)
    else:
        return variable


def get_first_available_time(variables: List[str]) -> float:
    first_otimes = []

    for variable in variables:
        files = glob.glob('{0}.[0-9]*'.format(get_source_variable(variable)))
        if len(files) == 0:
            raise Exception('Cannot find any average files for {0}'.format(variable))

        otimes = [int(Path(f).name.split('.')[-1]) for f in files]
        first_otimes.append(min(otimes))

    return max(first_otimes) * 10**iotimeprec


@dataclass(frozen=True)
class DumpSeries:
    otimes: List[int]
    paths: Dict[int, Path]


@dataclass(frozen=True)
class OutputFile:
    handle: object
    time_var: object
    var: object
    start_index: int
    last_otime: Optional[int]

    def close(self) -> None:
        self.handle.close()


@dataclass(frozen=True)
class DumpOutputSpec:
    filename: str
    variable: str
    dim: collections.OrderedDict


@dataclass(frozen=True)
class ExistingOutputState:
    stored_otimes: List[Optional[int]]
    last_good_otime: Optional[int]


JOB_DIRS: Optional[List[Path]] = None


def _scratch_case_dir(casename: str) -> Path:
    user = os.environ.get('USER')
    if not user:
        raise RuntimeError('Environment variable USER is not set')
    return Path('/scratch/general/vast') / user / casename


def _parse_jobids(jobids: List[str]) -> List[str]:
    if not jobids:
        raise ValueError('--jobids requires at least one jobid')
    for jid in jobids:
        if not str(jid).isdigit():
            raise ValueError(f'Invalid jobid "{jid}" (expected integer)')
    return [str(j) for j in jobids]


def _parse_dump_filename(filename: str, variable: str) -> Optional[int]:
    base = Path(filename).name
    parts = base.split('.')
    if len(parts) != 2 or parts[0] != variable:
        return None
    try:
        return int(parts[1])
    except ValueError:
        return None


def _discover_dump_series(variable: str, job_dirs: List[Path]) -> DumpSeries:
    paths: Dict[int, Path] = {}
    pattern = f'{glob.escape(variable)}.*'

    for job_dir in job_dirs:
        for f in glob.glob(str(job_dir / pattern)):
            otime = _parse_dump_filename(f, variable)
            if otime is None:
                continue
            paths[otime] = Path(f)

    if not paths:
        raise Exception(f'Cannot find any 3D dump files for "{variable}"')

    return DumpSeries(otimes=sorted(paths.keys()), paths=paths)


def _build_local_otimes() -> List[int]:
    otimes: List[int] = []
    for t, time in enumerate(np.arange(starttime, endtime + sampletime, sampletime)):
        otime = int(round(time / 10**iotimeprec))
        if doubledump and t > 0:
            with open(f'time.{otime:07d}', 'rb') as infile:
                timedata = struct.unpack('=QQi', infile.read())
            otime2 = int(round((timedata[0] - timedata[1]) * 10**(-iotimeprec - 9) - 0.5))
            otimes.append(otime2)
        otimes.append(otime)
    return otimes


def _build_output_dims(variable: str) -> collections.OrderedDict:
    dim = collections.OrderedDict()
    dim['time'] = []
    dim['z'] = range(kmax)
    dim['y'] = range(jtot)
    dim['x'] = range(itot)

    base_variable = variable[:-4] if variable.endswith('_avg') else variable
    if base_variable == 'u':
        dim['xh'] = dim.pop('x')
    if base_variable == 'v':
        dim['yh'] = dim.pop('y')
    if base_variable in HALF_LEVEL_VARS:
        dim['zh'] = dim.pop('z')

    return dim


def _sorted_dims(dim_keys) -> tuple:
    ordered_dims = ['time', 'z', 'zh', 'y', 'yh', 'x', 'xh']
    return tuple(value for value in ordered_dims if value in dim_keys)


def _expected_dim_sizes(dim) -> Dict[str, int]:
    sizes = {}
    for key, value in dim.items():
        if key == 'time':
            continue
        sizes[key] = len(value)
    return sizes


def _validate_existing_output(ncfile, filename: str, variable: str, dim):
    if 'time' not in ncfile.variables:
        raise Exception(f'Existing file "{filename}" is missing the "time" variable')
    if variable not in ncfile.variables:
        raise Exception(f'Existing file "{filename}" is missing the "{variable}" variable')

    expected_dims = _sorted_dims(dim.keys())
    actual_dims = tuple(ncfile.variables[variable].dimensions)
    if actual_dims != expected_dims:
        raise Exception(
            f'Existing file "{filename}" has dimensions {actual_dims}, expected {expected_dims}'
        )

    for key, expected_size in _expected_dim_sizes(dim).items():
        if key not in ncfile.dimensions:
            raise Exception(f'Existing file "{filename}" is missing dimension "{key}"')
        actual_size = len(ncfile.dimensions[key])
        if actual_size != expected_size:
            raise Exception(
                f'Existing file "{filename}" has dimension "{key}" size {actual_size}, '
                f'expected {expected_size}'
            )

    return ncfile.variables['time'], ncfile.variables[variable]


def _open_existing_output(filename: str, variable: str, dim) -> OutputFile:
    ncfile = mht.nc.Dataset(filename, 'a')
    try:
        time_var, var = _validate_existing_output(ncfile, filename, variable, dim)
        start_index = 0
        last_otime = None
        for idx in range(len(time_var) - 1, -1, -1):
            if _is_record_complete(time_var, var, idx):
                start_index = idx + 1
                last_time = float(time_var[idx])
                last_otime = int(round(last_time / 10**iotimeprec))
                break

        return OutputFile(
            handle=ncfile,
            time_var=time_var,
            var=var,
            start_index=start_index,
            last_otime=last_otime,
        )
    except Exception:
        ncfile.close()
        raise


def _create_new_output(filename: str, variable: str, dim) -> OutputFile:
    ncfile = mht.Create_ncfile(grid, filename, variable, dim, precision, compression)
    return OutputFile(
        handle=ncfile,
        time_var=ncfile.dimvar['time'],
        var=ncfile.var,
        start_index=0,
        last_otime=None,
    )


def _is_record_complete(time_var, var, index: int) -> bool:
    time_value = time_var[index]
    if np.ma.is_masked(time_value):
        return False
    data = var[index]
    if np.ma.isMaskedArray(data):
        return not bool(np.ma.getmaskarray(data).any())
    return True


def _inspect_existing_output(filename: str, variable: str, dim) -> ExistingOutputState:
    ncfile = mht.nc.Dataset(filename, 'r')
    try:
        time_var, var = _validate_existing_output(ncfile, filename, variable, dim)
        stored_otimes: List[Optional[int]] = []
        for time_value in time_var[:]:
            if np.ma.is_masked(time_value):
                stored_otimes.append(None)
            else:
                stored_otimes.append(int(round(float(time_value) / 10**iotimeprec)))
        last_good_otime = None
        for idx in range(len(stored_otimes) - 1, -1, -1):
            if _is_record_complete(time_var, var, idx):
                last_good_otime = stored_otimes[idx]
                break
        return ExistingOutputState(
            stored_otimes=stored_otimes,
            last_good_otime=last_good_otime,
        )
    finally:
        ncfile.close()


def _truncate_existing_output(spec: DumpOutputSpec, keep_count: int) -> None:
    tmp_filename = f'{spec.filename}.truncate.tmp'
    if os.path.exists(tmp_filename):
        os.remove(tmp_filename)

    src = mht.nc.Dataset(spec.filename, 'r')
    output = None
    try:
        time_var, var = _validate_existing_output(src, spec.filename, spec.variable, spec.dim)
        output = _create_new_output(tmp_filename, spec.variable, spec.dim)
        if keep_count > 0:
            output.time_var[:keep_count] = time_var[:keep_count]
            output.var[:keep_count] = var[:keep_count]
        output.close()
        output = None
        src.close()
        src = None
        os.replace(tmp_filename, spec.filename)
    finally:
        if output is not None:
            output.close()
        if src is not None:
            src.close()
        if os.path.exists(tmp_filename):
            os.remove(tmp_filename)


def _truncate_outputs_if_requested() -> None:
    specs = [
        DumpOutputSpec(
            filename=f'{get_source_variable(variable)}.nc',
            variable=get_source_variable(variable),
            dim=_build_output_dims(get_source_variable(variable)),
        )
        for variable in variables
    ]

    existing_specs = [spec for spec in specs if os.path.isfile(spec.filename)]
    if not existing_specs:
        print('No existing 3D netCDF files found to truncate.')
        return

    states = {
        spec.filename: _inspect_existing_output(spec.filename, spec.variable, spec.dim)
        for spec in existing_specs
    }

    if any(state.last_good_otime is None for state in states.values()):
        common_last_good_otime = None
    else:
        common_last_good_otime = min(
            state.last_good_otime for state in states.values() if state.last_good_otime is not None)

    if common_last_good_otime is None:
        print('Truncating 3D outputs to zero timesteps (no common complete timestep found).')
    else:
        print(f'Truncating 3D outputs to common last complete otime {common_last_good_otime:07d}.')

    for spec in existing_specs:
        state = states[spec.filename]
        if common_last_good_otime is None:
            keep_count = 0
        else:
            keep_count = sum(
                1
                for otime in state.stored_otimes
                if otime is not None and otime <= common_last_good_otime
            )
        _truncate_existing_output(spec, keep_count)


def convert_to_nc(variables):
    for variable in variables:
        output = None
        try:
            source_variable = get_source_variable(variable)
            filename = f'{source_variable}.nc'
            dim = _build_output_dims(source_variable)

            if JOB_DIRS is None:
                otimes = _build_local_otimes()
                series = None
            else:
                series = _discover_dump_series(source_variable, JOB_DIRS)
                otimes = list(series.otimes)

            if append and os.path.isfile(filename):
                output = _open_existing_output(filename, source_variable, dim)
                if output.last_otime is not None:
                    otimes = [otime for otime in otimes if otime > output.last_otime]
                if not otimes:
                    print(f'No new timesteps for {filename}; skipping append.')
                    output.close()
                    output = None
                    continue
                time_offset = output.start_index
            else:
                output = _create_new_output(filename, source_variable, dim)
                time_offset = 0

            for tout, otime in enumerate(otimes, start=time_offset):
                if JOB_DIRS is None:
                    f_path = f'{source_variable}.{otime:07d}'
                else:
                    try:
                        f_path = str(series.paths[otime])  # type: ignore[union-attr]
                    except KeyError:
                        raise Exception(f'Missing 3D dump file for "{source_variable}" time {otime:07d}')

                try:
                    fin = mht.Read_binary(grid, f_path)
                except Exception as ex:
                    raise Exception(f'Cannot read 3D dump file "{f_path}": {ex}')

                print("Processing %8s, time=%7i" % (source_variable, otime))
                if perslice:
                    for k in range(kmax):
                        output.var[tout, k, :, :] = fin.read(itot * jtot)
                else:
                    output.var[tout, :, :, :] = fin.read(itot * jtot * kmax)
                output.time_var[tout] = otime * 10**iotimeprec
                fin.close()

            output.close()
            output = None
        except Exception as ex:
            if output is not None:
                output.close()
            print(ex)
            print("Failed to create %s" % filename)


parser = argparse.ArgumentParser(
    description='Convert MicroHH 3D binary to netCDF4 files.')
parser.add_argument('-d', '--directory', help='directory')
parser.add_argument('-f', '--filename', help='ini file name')
parser.add_argument('-v', '--vars', nargs='*', help='variable names')
parser.add_argument(
    '-p',
    '--precision',
    help='precision',
    choices=[
        'single',
        'double'])
parser.add_argument(
    '-o',
    '--order',
    help='order',
    choices=[
        2, 4], type=int)
parser.add_argument(
    '-t0',
    '--starttime',
    help='first time step to be parsed',
    type=float)
parser.add_argument(
    '-t1',
    '--endtime',
    help='last time step to be parsed',
    type=float)
parser.add_argument(
    '-tstep',
    '--sampletime',
    help='time interval to be parsed',
    type=float)
parser.add_argument(
    '-s',
    '--perslice',
    help='read/write per horizontal slice',
    action='store_true')
parser.add_argument(
    '-c',
    '--nocompression',
    help='do not compress the netcdf file',
    action='store_true')
parser.add_argument(
    '-kmax',
    '--kmax',
    help='reduce vertical extent 3D files',
    type=int)
parser.add_argument('-n', '--nprocs', help='Number of processes', type=int)
parser.add_argument(
    '-cn',
    '--casename',
    help='Case name under /scratch/general/vast/$USER/<casename>/')
parser.add_argument(
    '-j',
    '--jobids',
    nargs='+',
    help='Restart chain jobids under /scratch/general/vast/$USER/<casename>/<jobid>/')
parser.add_argument(
    '-a',
    '--append',
    help='append only new timesteps to existing netCDF files when possible',
    action='store_true')
parser.add_argument(
    '--truncate',
    help='truncate existing netCDF files back to the common last complete timestep before optional append',
    action='store_true')
parser.add_argument(
    '--avg',
    help='target running-average 3D fields',
    action='store_true')

args = parser.parse_args()

if args.jobids is not None:
    if args.casename is None:
        print('ERROR: --casename is required when using --jobids', file=sys.stderr)
        sys.exit(2)

    jobids = _parse_jobids(args.jobids)
    case_dir = _scratch_case_dir(args.casename)
    JOB_DIRS = [case_dir / jid for jid in jobids]

    for d in JOB_DIRS:
        if not d.is_dir():
            print(f'ERROR: Job directory not found: {d}', file=sys.stderr)
            sys.exit(2)

    out_dir = JOB_DIRS[-1]
    if args.directory is not None and Path(args.directory).resolve() != out_dir.resolve():
        print(f'WARNING: Ignoring --directory "{args.directory}" because --jobids is used.')
    os.chdir(out_dir)

    if args.starttime is not None or args.endtime is not None or args.sampletime is not None:
        print('WARNING: Ignoring --starttime/--endtime/--sampletime because --jobids is used.')
else:
    if args.directory is not None:
        os.chdir(args.directory)

nl = mht.Read_namelist(args.filename)
itot = nl['grid']['itot']
jtot = nl['grid']['jtot']
ktot = nl['grid']['ktot']
kmax = args.kmax if args.kmax is not None else ktot
kmax = min(kmax, ktot)

starttime = args.starttime if args.starttime is not None else nl['time']['starttime']
endtime = args.endtime if args.endtime is not None else nl['time']['endtime']
avg = args.avg
if avg:
    starttime = args.starttime if args.starttime is not None else nl['average'].get('starttime', 0.)
    sampletime = args.sampletime if args.sampletime is not None else nl['average']['sampletime']
else:
    sampletime = args.sampletime if args.sampletime is not None else nl['dump']['sampletime']
try:
    doubledump = (nl['dump']['swdoubledump'] == 1) and not avg
except Exception:
    doubledump = False

try:
    iotimeprec = nl['time']['iotimeprec']
except KeyError:
    iotimeprec = 0.

variables = args.vars if args.vars is not None else (
        nl['average']['averagelist'] if avg else nl['dump']['dumplist'])
if isinstance(variables, str):
    variables = [variables]

if avg and JOB_DIRS is None and args.starttime is None:
    starttime = get_first_available_time(variables)

precision = args.precision
perslice = args.perslice
compression = not(args.nocompression)
append = args.append
truncate = args.truncate
nprocs = args.nprocs if args.nprocs is not None else len(variables)

try:
    order = args.order if args.order is not None else nl['grid']['swspatialorder']
except KeyError:
    order = 2

if JOB_DIRS is None:
    for time in np.arange(starttime, endtime, sampletime):
        otime = int(round(time / 10**iotimeprec))
        source_variables = [get_source_variable(variable) for variable in variables]
        if not any(glob.glob('{0}.{1:07d}'.format(variable, otime)) for variable in source_variables):
            endtime = time - sampletime
            break

grid = mht.Read_grid(itot, jtot, ktot, order=order)

if kmax < ktot:
    grid.dim['z'] = grid.dim['z'][:kmax]
    grid.dim['zh'] = grid.dim['zh'][:kmax+1]

if truncate:
    _truncate_outputs_if_requested()
    if not append:
        sys.exit(0)

chunks = [variables[i::nprocs] for i in range(nprocs)]

pool = Pool(processes=nprocs)
for _ in pool.imap_unordered(convert_to_nc, chunks):
    pass

pool.close()
pool.join()
