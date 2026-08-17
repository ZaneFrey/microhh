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
import sys
import microhh_tools as mht  # available in microhh/python directory
import argparse
import collections
import glob
import numpy as np
from multiprocessing import Pool, set_start_method
import platform
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

if platform.system() == 'Darwin':
    set_start_method('fork')

HALF_LEVEL_VARS = {'w', 'lflx', 'sflx'}


@dataclass(frozen=True)
class CrossSeries:
    at_surface: bool
    halflevel: str
    indices: List[int]
    otimes: List[int]
    paths: Dict[Tuple[Optional[int], int], Path]


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
class CrossOutputSpec:
    filename: str
    variable: str
    mode: str
    dim: collections.OrderedDict
    indexes_local: List[int]


@dataclass(frozen=True)
class ExistingOutputState:
    stored_otimes: List[Optional[int]]
    last_good_otime: Optional[int]


JOB_DIRS: Optional[List[Path]] = None


def get_source_variable(variable: str) -> str:
    if avg and not variable.endswith('_avg'):
        return '{0}_avg'.format(variable)
    else:
        return variable


def get_base_variable(variable: str) -> str:
    if variable.endswith('_avg'):
        return variable[:-4]
    else:
        return variable


def get_halflevel(variable: str) -> str:
    base_variable = get_base_variable(variable)
    return ''.join([
        '1' if base_variable == 'u' else '0',
        '1' if base_variable == 'v' else '0',
        '1' if base_variable in HALF_LEVEL_VARS else '0'])


def get_vertical_coordinate(variable: str):
    if get_halflevel(variable)[2] == '1':
        return grid.dim['zh'][:-1]
    else:
        return grid.dim['z']


def get_cross_coordinate(mode: str, variable: str):
    halflevel = get_halflevel(variable)

    if mode == 'xy':
        return get_vertical_coordinate(variable)
    elif mode == 'xz':
        return grid.dim['yh'] if halflevel[1] == '1' else grid.dim['y']
    elif mode == 'yz':
        return grid.dim['xh'] if halflevel[0] == '1' else grid.dim['x']
    else:
        raise ValueError(f'Unsupported mode "{mode}"')


def get_average_indexes(mode: str, variable: str) -> List[int]:
    if indexes is not None:
        return indexes

    values = nl['cross'][mode]
    values = [values] if not isinstance(values, list) else values
    coordinate = get_cross_coordinate(mode, variable)

    return [int(np.abs(coordinate - value).argmin()) for value in values]


def get_first_available_time(variables: List[str]) -> float:
    first_otimes = []

    for variable in variables:
        files = glob.glob('{0}.[0-9]*'.format(get_source_variable(variable)))
        if len(files) == 0:
            raise Exception('Cannot find any average files for {0}'.format(variable))

        otimes = [int(Path(f).name.split('.')[-1]) for f in files]
        first_otimes.append(min(otimes))

    return max(first_otimes) * 10**iotimeprec


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


def _parse_cross_filename(
    filename: str, variable: str, mode: str
) -> Optional[Tuple[str, Optional[int], int]]:
    """
    Parse MicroHH cross-section binary filenames.

    Returns (halflevel, index|None, otime) or None if it does not match.
    """
    base = Path(filename).name
    parts = base.split('.')
    if len(parts) < 4:
        return None

    if parts[0] != variable or parts[1] != mode:
        return None

    # Surface: var.xy.<halflevel>.<otime>
    if len(parts) == 4:
        halflevel = parts[2]
        try:
            otime = int(parts[3])
        except ValueError:
            return None
        return halflevel, None, otime

    # Volume: var.<mode>.<halflevel>.<index>.<otime>
    if len(parts) == 5:
        halflevel = parts[2]
        try:
            index = int(parts[3])
            otime = int(parts[4])
        except ValueError:
            return None
        return halflevel, index, otime

    return None


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


def _discover_cross_series(
    variable: str,
    mode: str,
    job_dirs: List[Path],
    indexes_requested: Optional[List[int]],
    iotimeprec: float,
) -> CrossSeries:
    """
    Discover available cross-section binaries across multiple job directories.

    - No copying/moving: only returns paths to read.
    - If duplicates exist, later job_dirs override earlier ones.
    """
    # Build a key->path map, overwriting with later jobids to prefer the latest.
    paths: Dict[Tuple[Optional[int], int], Path] = {}
    halflevels_surface: Set[str] = set()
    halflevels_volume: Set[str] = set()
    has_surface = False
    has_volume = False

    pattern = f"{glob.escape(variable)}.{glob.escape(mode)}.*"
    for job_dir in job_dirs:
        for f in glob.glob(str(job_dir / pattern)):
            parsed = _parse_cross_filename(f, variable, mode)
            if parsed is None:
                continue
            halflevel, index, otime = parsed

            # Only 'xy' can have a surface file (no index).
            if index is None and mode != 'xy':
                continue

            key = (index, otime)
            paths[key] = Path(f)

            if index is None:
                has_surface = True
                halflevels_surface.add(halflevel)
            else:
                has_volume = True
                halflevels_volume.add(halflevel)

    if not paths:
        raise Exception(f'Cannot find any cross-section files for "{variable}" mode "{mode}"')

    at_surface = (mode == 'xy') and has_surface
    if at_surface and has_volume:
        # Match the original script's behavior: if a surface file exists, treat it as surface output.
        print(
            f'WARNING: Found both surface and volume cross files for "{variable}" (mode "{mode}"); '
            f'using surface files.'
        )

    if at_surface:
        if len(halflevels_surface) != 1:
            raise Exception(
                f'Ambiguous halflevel for surface files of "{variable}" (mode "{mode}"): '
                f'{sorted(halflevels_surface)}'
            )
        halflevel = next(iter(halflevels_surface))
        indices_local = [-1]
        otimes = sorted({otime for (idx, otime) in paths.keys() if idx is None})
        # Restrict to surface-only keys.
        paths = {(idx, otime): p for (idx, otime), p in paths.items() if idx is None}
        return CrossSeries(
            at_surface=True, halflevel=halflevel, indices=indices_local, otimes=otimes, paths=paths
        )

    # Volume files:
    if len(halflevels_volume) != 1:
        raise Exception(
            f'Ambiguous halflevel for volume files of "{variable}" (mode "{mode}"): '
            f'{sorted(halflevels_volume)}'
        )
    halflevel = next(iter(halflevels_volume))

    if indexes_requested is not None:
        indices_local = list(indexes_requested)
    else:
        indices_local = sorted({idx for (idx, _otime) in paths.keys() if idx is not None})
        if not indices_local:
            raise Exception(f'Cannot find any indexed cross-section files for "{variable}" mode "{mode}"')

    # Only keep volume keys.
    paths = {(idx, otime): p for (idx, otime), p in paths.items() if idx is not None}

    # Find otimes where all indices exist (skip incomplete times).
    time_sets: List[Set[int]] = []
    for idx in indices_local:
        tset = {otime for (i, otime) in paths.keys() if i == idx}
        if not tset:
            raise Exception(f'Cannot find any times for "{variable}" mode "{mode}" index {idx}')
        time_sets.append(tset)

    otimes = sorted(set.intersection(*time_sets)) if time_sets else []
    if not otimes:
        raise Exception(
            f'No complete time steps found for "{variable}" mode "{mode}" across indices {indices_local}'
        )

    # Sanity: iotimeprec should be usable even though we discover from files
    _ = otimes[0] * (10 ** iotimeprec)

    return CrossSeries(
        at_surface=False, halflevel=halflevel, indices=indices_local, otimes=otimes, paths=paths
    )


def _build_output_dims(
    at_surface: bool,
    mode: str,
    halflevel: str,
    indexes_local: List[int],
) -> Tuple[collections.OrderedDict, int]:
    dim = collections.OrderedDict()
    dim['time'] = []
    dim['z'] = range(ktot)
    dim['y'] = range(jtot)
    dim['x'] = range(itot)

    if at_surface:
        dim.pop('z')
        n = itot * jtot
    elif mode == 'xy':
        dim.update({'z': []})
        n = itot * jtot
    elif mode == 'xz':
        dim.update({'y': []})
        n = itot * ktot
    elif mode == 'yz':
        dim.update({'x': []})
        n = ktot * jtot
    else:
        raise ValueError(f'Unsupported mode "{mode}"')

    if halflevel[0] == '1':
        dim['xh'] = dim.pop('x')
    if halflevel[1] == '1':
        dim['yh'] = dim.pop('y')
    if halflevel[2] == '1':
        dim['zh'] = dim.pop('z')

    return dim, n


def _sorted_dims(dim_keys) -> Tuple[str, ...]:
    ordered_dims = ['time', 'z', 'zh', 'y', 'yh', 'x', 'xh']
    return tuple(value for value in ordered_dims if value in dim_keys)


def _expected_dim_sizes(dim, indexes_local: List[int]) -> Dict[str, int]:
    sizes = {}
    for key, value in dim.items():
        if key == 'time':
            continue
        sizes[key] = len(indexes_local) if value == [] else len(value)
    return sizes


def _validate_existing_output(ncfile, filename: str, variable: str, dim, indexes_local: List[int]):
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

    for key, expected_size in _expected_dim_sizes(dim, indexes_local).items():
        if key not in ncfile.dimensions:
            raise Exception(f'Existing file "{filename}" is missing dimension "{key}"')
        actual_size = len(ncfile.dimensions[key])
        if actual_size != expected_size:
            raise Exception(
                f'Existing file "{filename}" has dimension "{key}" size {actual_size}, '
                f'expected {expected_size}'
            )

    return ncfile.variables['time'], ncfile.variables[variable]


def _open_existing_output(
    filename: str,
    variable: str,
    dim,
    indexes_local: List[int],
    iotimeprec: float,
) -> OutputFile:
    ncfile = mht.nc.Dataset(filename, 'a')
    try:
        time_var, var = _validate_existing_output(
            ncfile, filename, variable, dim, indexes_local)
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


def _create_new_output(
    filename: str,
    variable: str,
    dim,
    indexes_local: List[int],
) -> OutputFile:
    ncfile = mht.Create_ncfile(
        grid, filename, variable, dim, precision, compression)

    for key, val in dim.items():
        if key == 'time':
            continue
        elif val == []:
            ncfile.dimvar[key][:] = grid.dim[key][indexes_local]

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


def _inspect_existing_output(
    filename: str,
    variable: str,
    dim,
    indexes_local: List[int],
) -> ExistingOutputState:
    ncfile = mht.nc.Dataset(filename, 'r')
    try:
        time_var, var = _validate_existing_output(
            ncfile, filename, variable, dim, indexes_local)
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


def _truncate_existing_output(
    spec: CrossOutputSpec,
    keep_count: int,
) -> None:
    tmp_filename = f'{spec.filename}.truncate.tmp'
    if os.path.exists(tmp_filename):
        os.remove(tmp_filename)

    src = mht.nc.Dataset(spec.filename, 'r')
    output = None
    try:
        time_var, var = _validate_existing_output(
            src, spec.filename, spec.variable, spec.dim, spec.indexes_local)
        output = _create_new_output(
            filename=tmp_filename,
            variable=spec.variable,
            dim=spec.dim,
            indexes_local=spec.indexes_local,
        )
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


def _resolve_output_spec(variable: str, mode: str) -> Optional[CrossOutputSpec]:
    source_variable = get_source_variable(variable)

    if avg:
        at_surface = False
        filename = "{0}.{1}.nc".format(source_variable, mode)
        halflevel = get_halflevel(source_variable)
        indexes_local = get_average_indexes(mode, source_variable)
        dim, _ = _build_output_dims(at_surface, mode, halflevel, indexes_local)
        return CrossOutputSpec(
            filename=filename,
            variable=source_variable,
            mode=mode,
            dim=dim,
            indexes_local=indexes_local,
        )

    if JOB_DIRS is None:
        otime = int(round(starttime / 10**iotimeprec))

        if os.path.isfile("{0}.xy.000.{1:07d}".format(source_variable, otime)):
            if mode != 'xy':
                return None
            at_surface = True
        else:
            at_surface = False

        filename = "{0}.{1}.nc".format(variable, mode)
        halflevel = '000'
        if at_surface:
            indexes_local = [-1]
        else:
            if indexes is None:
                indexes_local, halflevel = mht.get_cross_indices(source_variable, mode)
            else:
                indexes_local = indexes
                files = glob.glob("{0:}.{1}.*.{2:05d}.{3:07d}".format(
                        source_variable, mode, indexes_local[0], otime))
                if len(files) == 0:
                    raise Exception('Cannot find any cross-section')
                halflevel = files[0].split('.')[-3]
    else:
        series = _discover_cross_series(
            variable=source_variable,
            mode=mode,
            job_dirs=JOB_DIRS,
            indexes_requested=indexes,
            iotimeprec=iotimeprec,
        )
        at_surface = series.at_surface
        if at_surface and mode != 'xy':
            return None
        filename = "{0}.{1}.nc".format(source_variable, mode)
        halflevel = series.halflevel
        indexes_local = series.indices

    dim, _ = _build_output_dims(at_surface, mode, halflevel, indexes_local)
    return CrossOutputSpec(
        filename=filename,
        variable=source_variable,
        mode=mode,
        dim=dim,
        indexes_local=indexes_local,
    )


def _truncate_outputs_if_requested() -> None:
    specs: List[CrossOutputSpec] = []
    for variable in variables:
        for mode in modes:
            spec = _resolve_output_spec(variable, mode)
            if spec is not None:
                specs.append(spec)

    existing_specs = [spec for spec in specs if os.path.isfile(spec.filename)]
    if not existing_specs:
        print('No existing cross netCDF files found to truncate.')
        return

    states = {
        spec.filename: _inspect_existing_output(
            spec.filename, spec.variable, spec.dim, spec.indexes_local)
        for spec in existing_specs
    }

    if any(state.last_good_otime is None for state in states.values()):
        common_last_good_otime = None
    else:
        common_last_good_otime = min(
            state.last_good_otime for state in states.values() if state.last_good_otime is not None)

    if common_last_good_otime is None:
        print('Truncating cross outputs to zero timesteps (no common complete timestep found).')
    else:
        print(f'Truncating cross outputs to common last complete otime {common_last_good_otime:07d}.')

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
    # Loop over the different variables and crosssections
    for variable in variables:
        source_variable = get_source_variable(variable)
        for mode in modes:
            output = None
            try:
                if avg:
                    at_surface = False
                    filename = "{0}.{1}.nc".format(source_variable, mode)
                    halflevel = get_halflevel(source_variable)
                    indexes_local = get_average_indexes(mode, source_variable)

                    if JOB_DIRS is None:
                        times = list(np.arange(starttime, endtime + sampletime, sampletime))
                        series = None
                    else:
                        series = _discover_dump_series(source_variable, JOB_DIRS)
                        times = [ot * (10 ** iotimeprec) for ot in series.otimes]
                else:
                    if JOB_DIRS is None:
                        otime = int(round(starttime / 10**iotimeprec))

                        if os.path.isfile("{0}.xy.000.{1:07d}".format(source_variable, otime)):
                            if mode != 'xy':
                                continue
                            at_surface = True
                        else:
                            at_surface = False

                        filename = "{0}.{1}.nc".format(source_variable, mode)
                        halflevel = '000'
                        if not at_surface:
                            if indexes is None:
                                indexes_local, halflevel = mht.get_cross_indices(source_variable, mode)
                            else:
                                indexes_local = indexes

                                files = glob.glob("{0:}.{1}.*.{2:05d}.{3:07d}".format(
                                        source_variable, mode, indexes_local[0], otime))
                                if len(files) == 0:
                                    raise Exception('Cannot find any cross-section')
                                halflevel = files[0].split('.')[-3]

                        times = list(np.arange(starttime, endtime + sampletime, sampletime))
                        series = None
                    else:
                        series = _discover_cross_series(
                            variable=source_variable,
                            mode=mode,
                            job_dirs=JOB_DIRS,
                            indexes_requested=indexes,
                            iotimeprec=iotimeprec,
                        )
                        at_surface = series.at_surface
                        if at_surface and mode != 'xy':
                            continue
                        filename = "{0}.{1}.nc".format(source_variable, mode)
                        halflevel = series.halflevel
                        indexes_local = series.indices
                        times = [ot * (10 ** iotimeprec) for ot in series.otimes]

                dim, n = _build_output_dims(at_surface, mode, halflevel, indexes_local)
                time_offset = 0

                if append and os.path.isfile(filename):
                    output = _open_existing_output(
                        filename=filename,
                        variable=source_variable,
                        dim=dim,
                        indexes_local=indexes_local,
                        iotimeprec=iotimeprec,
                    )
                    time_offset = output.start_index
                    if output.last_otime is not None:
                        times = [
                            time
                            for time in times
                            if int(round(time / 10**iotimeprec)) > output.last_otime
                        ]
                    if not times:
                        print(f'No new timesteps for {filename}; skipping append.')
                        output.close()
                        output = None
                        continue
                else:
                    output = _create_new_output(
                        filename=filename,
                        variable=source_variable,
                        dim=dim,
                        indexes_local=indexes_local,
                    )

                for t, time in enumerate(times, start=time_offset):
                    if avg:
                        otime = int(round((time) / 10**iotimeprec))

                        if JOB_DIRS is None:
                            f_path = "{0:}.{1:07d}".format(source_variable, otime)
                        else:
                            try:
                                f_path = str(series.paths[otime])  # type: ignore[union-attr]
                            except KeyError:
                                raise Exception(
                                    f'Missing 3D average file for "{source_variable}" time {otime:07d}'
                                )

                        try:
                            fin = mht.Read_binary(grid, str(f_path))
                        except Exception as ex:
                            print(ex)
                            break

                        print("Processing %8s, time=%7i" % (source_variable, otime))

                        field = fin.read(itot * jtot * ktot).reshape((ktot, jtot, itot))
                        time_block = np.empty(_time_block_shape(mode, len(indexes_local), at_surface), dtype=grid.dtype)

                        for k, index in enumerate(indexes_local):
                            if mode == 'xy':
                                time_block[k, :, :] = field[index, :, :]
                            elif mode == 'xz':
                                time_block[:, k, :] = field[:, index, :]
                            elif mode == 'yz':
                                time_block[:, :, k] = field[:, :, index]

                        output.var[t, ...] = time_block
                        output.time_var[t] = time
                        fin.close()
                        continue

                    time_complete = True
                    for k in range(len(indexes_local)):
                        index = indexes_local[k]
                        otime = int(round((time) / 10**iotimeprec))

                        if JOB_DIRS is None:
                            if at_surface:
                                f_in = "{0}.{1}.{2}.{3:07d}".format(
                                    source_variable, mode, halflevel, otime)
                            else:
                                f_in = "{0:}.{1}.{2}.{3:05d}.{4:07d}".format(
                                    source_variable, mode, halflevel, index, otime)
                            f_path = f_in
                        else:
                            if at_surface:
                                key = (None, otime)
                            else:
                                key = (index, otime)
                            try:
                                f_path = series.paths[key]  # type: ignore[union-attr]
                            except KeyError:
                                raise Exception(
                                    f'Missing cross file for "{source_variable}" mode "{mode}" '
                                    f'index {index} time {otime:07d}'
                                )
                        try:
                            fin = mht.Read_binary(grid, str(f_path))
                        except Exception as ex:
                            print(ex)
                            time_complete = False
                            break

                        print(
                            "Processing %8s, time=%7i, index=%4i" %
                            (source_variable, otime, index))

                        if at_surface:
                            output.var[t, :, :] = fin.read(n)
                        elif mode == 'xy':
                            output.var[t, k, :, :] = fin.read(n)
                        elif mode == 'xz':
                            output.var[t, :, k, :] = fin.read(n)
                        elif mode == 'yz':
                            output.var[t, :, :, k] = fin.read(n)

                        fin.close()
                    if not time_complete:
                        break
                    output.time_var[t] = time
                output.close()
                output = None

            except Exception as ex:
                if output is not None:
                    output.close()
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
    modes = [m for m in cross_modes if m in nl['cross']]

    # Check if there are paths in the cross-list
    if 'xy' not in modes:
        for v in nl['cross']['crosslist']:
            if isinstance(v, dict) and 'path' in v:
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

if avg and JOB_DIRS is None and args.starttime is None:
    starttime = get_first_available_time(variables)

precision = args.precision
nprocs = args.nprocs if args.nprocs is not None else len(variables)
compression = not(args.nocompression)
append = args.append
truncate = args.truncate
try:
    order = args.order if args.order is not None else nl['grid']['swspatialorder']
except KeyError:
    order = 2

# End option parsing
grid = mht.Read_grid(itot, jtot, ktot, order = order)

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
