#!/usr/bin/env python3
"""Reproduce Calaf et al. (2010) diagnostics from the MicroHH case output.

The script makes paper-style plots of instantaneous velocity, the mean wind
profile, momentum stresses, mean-energy conversion and transport, and turbine
power.  It also writes profile CSV files and a text summary of the principal
nondimensional statistics.

By default, the latter half of the available run is used as the averaging
window.  This is a pragmatic choice for the short demonstration run, not a
claim of statistical stationarity.  Use --start-time and --end-time to select
the averaging interval explicitly for a production run.

Examples
--------
Run from the case directory using the latter half of the available data::

    python3 calaf2010_analysis.py

Select an explicit averaging window::

    python3 calaf2010_analysis.py --start-time 21600 --end-time 43200

Skip the relatively expensive 3-D dispersive-stress calculation::

    python3 calaf2010_analysis.py --skip-3d
"""

from __future__ import annotations

import argparse
from configparser import ConfigParser
from dataclasses import dataclass
import os
from pathlib import Path
import tempfile
import warnings

import netCDF4 as nc
import numpy as np


# Keep Matplotlib's cache out of the source tree on systems with a read-only
# home directory.  This must be set before importing matplotlib.
_mpl_cache = Path(tempfile.gettempdir()) / "calaf2010-matplotlib"
_mpl_cache.mkdir(parents=True, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(_mpl_cache))

import matplotlib.pyplot as plt  # noqa: E402


KAPPA = 0.4
SHOW_FIGURES = False


@dataclass(frozen=True)
class CaseConfig:
    case_dir: Path
    name: str
    nx: int
    ny: int
    nz: int
    lx: float
    ly: float
    height: float
    diameter: float
    hub_height: float
    roughness: float
    dpdx: float

    @property
    def ustar(self) -> float:
        """Pressure-gradient friction velocity used by Calaf et al."""

        return float(np.sqrt(-self.dpdx * self.height))

    @property
    def turbine_bottom(self) -> float:
        return self.hub_height - 0.5 * self.diameter

    @property
    def turbine_top(self) -> float:
        return self.hub_height + 0.5 * self.diameter


@dataclass
class ProfileData:
    time: np.ndarray
    z: np.ndarray
    zh: np.ndarray
    selected: np.ndarray
    weights: np.ndarray
    u: np.ndarray
    u2: np.ndarray
    resolved_flux: np.ndarray
    diffusive_flux: np.ndarray
    total_flux: np.ndarray
    ustar_surface: np.ndarray
    rhoref: np.ndarray


@dataclass
class StressData:
    resolved: np.ndarray
    diffusive: np.ndarray
    total: np.ndarray
    turbulent: np.ndarray
    dispersive: np.ndarray
    separation_available: bool
    sampled_resolved: np.ndarray | None = None
    snapshot_count: int = 0


@dataclass
class TurbineData:
    time: np.ndarray
    ids: np.ndarray
    filtered_velocity: np.ndarray
    thrust: np.ndarray
    power_density: np.ndarray
    selected: np.ndarray
    weights: np.ndarray
    rho_hub: float
    plan_area_per_turbine: float


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Create Calaf et al. (2010)-style diagnostics from MicroHH output."
    )
    parser.add_argument(
        "--case-dir",
        type=Path,
        default=Path(__file__).resolve().parent,
        help="Case directory (default: directory containing this script).",
    )
    parser.add_argument(
        "--start-time",
        type=float,
        default=None,
        help="Start of averaging window in seconds (default: halfway through run).",
    )
    parser.add_argument(
        "--end-time",
        type=float,
        default=None,
        help="End of averaging window in seconds (default: last available record).",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help="Plot/output directory (default: CASE_DIR/calaf2010_figures).",
    )
    parser.add_argument(
        "--skip-3d",
        action="store_true",
        help="Skip 3-D time means and the turbulent/dispersive decomposition.",
    )
    parser.add_argument(
        "--skip-instantaneous",
        action="store_true",
        help="Skip the three instantaneous cross-section panels.",
    )
    parser.add_argument("--dpi", type=int, default=180, help="PNG resolution.")
    parser.add_argument("--show", action="store_true", help="Show figures interactively.")
    return parser.parse_args()


def read_case_config(case_dir: Path) -> CaseConfig:
    ini_path = case_dir / "calaf2010.ini"
    parser = ConfigParser()
    if not parser.read(ini_path):
        raise FileNotFoundError(f"Could not read {ini_path}")

    return CaseConfig(
        case_dir=case_dir,
        name="calaf2010",
        nx=parser.getint("grid", "itot"),
        ny=parser.getint("grid", "jtot"),
        nz=parser.getint("grid", "ktot"),
        lx=parser.getfloat("grid", "xsize"),
        ly=parser.getfloat("grid", "ysize"),
        height=parser.getfloat("grid", "zsize"),
        diameter=parser.getfloat("windfarm", "diameter"),
        hub_height=parser.getfloat("windfarm", "hubheight"),
        roughness=parser.getfloat("boundary", "z0m"),
        dpdx=parser.getfloat("force", "dpdx"),
    )


def as_float_array(variable: nc.Variable) -> np.ndarray:
    values = variable[:]
    if np.ma.isMaskedArray(values):
        values = values.filled(np.nan)
    return np.asarray(values, dtype=np.float64)


def select_time_indices(
    time: np.ndarray, start_time: float | None, end_time: float | None
) -> tuple[np.ndarray, float, float]:
    available_start = float(time[0])
    available_end = float(time[-1])
    end = available_end if end_time is None else min(float(end_time), available_end)
    if start_time is None:
        start = available_start + 0.5 * (end - available_start)
    else:
        start = max(float(start_time), available_start)

    selected = np.flatnonzero((time >= start - 1.0e-8) & (time <= end + 1.0e-8))
    if selected.size < 2:
        raise ValueError(
            f"Averaging window [{start}, {end}] contains fewer than two records "
            f"in the available interval [{available_start}, {available_end}]"
        )
    return selected, float(time[selected[0]]), float(time[selected[-1]])


def trapezoid_weights(time: np.ndarray) -> np.ndarray:
    if time.size == 1:
        return np.ones(1)
    dt = np.diff(time)
    if np.any(dt <= 0.0):
        raise ValueError("Time coordinate must be strictly increasing")
    weights = np.empty(time.size, dtype=np.float64)
    weights[0] = 0.5 * dt[0]
    weights[-1] = 0.5 * dt[-1]
    weights[1:-1] = 0.5 * (dt[:-1] + dt[1:])
    return weights / weights.sum()


def weighted_mean(values: np.ndarray, weights: np.ndarray) -> np.ndarray:
    return np.tensordot(weights, values, axes=(0, 0))


def read_profile_data(
    config: CaseConfig, start_time: float | None, end_time: float | None
) -> tuple[ProfileData, float, float]:
    paths = sorted(config.case_dir.glob(f"{config.name}.default.*.nc"))
    if not paths:
        raise FileNotFoundError("No MicroHH default statistics NetCDF file found")
    if len(paths) > 1:
        warnings.warn(
            "Multiple statistics files found; this script currently analyzes "
            f"the first one: {paths[0].name}"
        )

    with nc.Dataset(paths[0]) as dataset:
        group = dataset.groups["default"]
        time = as_float_array(dataset.variables["time"])
        selected, actual_start, actual_end = select_time_indices(
            time, start_time, end_time
        )
        weights = trapezoid_weights(time[selected])

        data = ProfileData(
            time=time,
            z=as_float_array(dataset.variables["z"]),
            zh=as_float_array(dataset.variables["zh"]),
            selected=selected,
            weights=weights,
            u=as_float_array(group.variables["u"]),
            u2=as_float_array(group.variables["u_2"]),
            resolved_flux=as_float_array(group.variables["u_w"]),
            diffusive_flux=as_float_array(group.variables["u_diff"]),
            total_flux=as_float_array(group.variables["u_flux"]),
            ustar_surface=as_float_array(group.variables["ustar"]),
            rhoref=as_float_array(group.variables["rhoref"]),
        )
    return data, actual_start, actual_end


def time_mean_profiles(profile: ProfileData) -> dict[str, np.ndarray | float]:
    idx = profile.selected
    weights = profile.weights
    return {
        "u": weighted_mean(profile.u[idx], weights),
        "u2": weighted_mean(profile.u2[idx], weights),
        "resolved_flux": weighted_mean(profile.resolved_flux[idx], weights),
        "diffusive_flux": weighted_mean(profile.diffusive_flux[idx], weights),
        "total_flux": weighted_mean(profile.total_flux[idx], weights),
        "ustar_surface": float(weighted_mean(profile.ustar_surface[idx], weights)),
    }


def collocate_u_at_w(u: np.ndarray) -> np.ndarray:
    """Interpolate u(z,y,xh) to the w(zh,y,x) locations on a uniform grid."""

    # xh locations are cell faces and x locations are cell centers.
    u_at_x = 0.5 * (u + np.roll(u, -1, axis=2))
    u_at_w = np.empty_like(u_at_x)
    u_at_w[0] = u_at_x[0]
    u_at_w[1:] = 0.5 * (u_at_x[:-1] + u_at_x[1:])
    return u_at_w


def calculate_dispersive_stress(
    config: CaseConfig, start_time: float, end_time: float
) -> tuple[np.ndarray, np.ndarray, np.ndarray, int]:
    u_path = config.case_dir / "u.nc"
    w_path = config.case_dir / "w.nc"
    if not u_path.exists() or not w_path.exists():
        raise FileNotFoundError("u.nc and w.nc are required for dispersive stress")

    with nc.Dataset(u_path) as u_ds, nc.Dataset(w_path) as w_ds:
        u_time = as_float_array(u_ds.variables["time"])
        w_time = as_float_array(w_ds.variables["time"])
        if u_time.shape != w_time.shape or not np.allclose(u_time, w_time):
            raise ValueError("u.nc and w.nc have different time coordinates")

        selected, _, _ = select_time_indices(u_time, start_time, end_time)
        weights = trapezoid_weights(u_time[selected])
        u_var = u_ds.variables["u"]
        w_var = w_ds.variables["w"]
        u_mean = np.zeros(u_var.shape[1:], dtype=np.float64)
        w_mean = np.zeros(w_var.shape[1:], dtype=np.float64)
        sampled_uw = np.zeros(w_var.shape[1], dtype=np.float64)

        print(f"Averaging {selected.size} paired 3-D u/w snapshots ...")
        for count, (record, weight) in enumerate(zip(selected, weights), start=1):
            u_now = np.asarray(u_var[record], dtype=np.float64)
            w_now = np.asarray(w_var[record], dtype=np.float64)
            u_mean += weight * u_now
            w_mean += weight * w_now
            sampled_uw += weight * np.mean(collocate_u_at_w(u_now) * w_now, axis=(1, 2))
            if count == 1 or count % 10 == 0 or count == selected.size:
                print(f"  3-D record {count:3d}/{selected.size:3d}")

        u_at_w = collocate_u_at_w(u_mean)
        u_plane = np.mean(u_at_w, axis=(1, 2))
        w_plane = np.mean(w_mean, axis=(1, 2))
        covariance = np.mean(
            (u_at_w - u_plane[:, None, None])
            * (w_mean - w_plane[:, None, None]),
            axis=(1, 2),
        )
        dispersive_stress = -covariance
        sampled_resolved_stress = -sampled_uw
        zh = as_float_array(w_ds.variables["zh"])
    return zh, dispersive_stress, sampled_resolved_stress, selected.size


def assemble_stresses(
    profile: ProfileData,
    means: dict[str, np.ndarray | float],
    config: CaseConfig,
    start_time: float,
    end_time: float,
    skip_3d: bool,
) -> StressData:
    resolved = -np.asarray(means["resolved_flux"]).copy()
    diffusive = -np.asarray(means["diffusive_flux"]).copy()
    total = -np.asarray(means["total_flux"]).copy()

    # The last staggered flux value is a numerical boundary/ghost value in
    # this MicroHH output, not a physical sample.  Impermeability makes all
    # vertical momentum fluxes vanish at the rigid lid.
    resolved[-1] = 0.0
    diffusive[-1] = 0.0
    total[-1] = 0.0

    if skip_3d:
        warnings.warn(
            "Skipping the 3-D mean fields: turbulent and dispersive stresses "
            "will not be separated."
        )
        nan = np.full_like(resolved, np.nan)
        return StressData(resolved, diffusive, total, resolved.copy(), nan, False)

    try:
        zh_3d, dispersive_3d, sampled_3d, count = calculate_dispersive_stress(
            config, start_time, end_time
        )
    except (FileNotFoundError, ValueError) as error:
        warnings.warn(f"Could not calculate dispersive stress: {error}")
        nan = np.full_like(resolved, np.nan)
        return StressData(resolved, diffusive, total, resolved.copy(), nan, False)

    dispersive = np.interp(profile.zh, zh_3d, dispersive_3d)
    sampled = np.interp(profile.zh, zh_3d, sampled_3d)
    dispersive[profile.zh > zh_3d[-1]] = np.nan
    sampled[profile.zh > zh_3d[-1]] = np.nan
    turbulent = resolved - dispersive
    return StressData(
        resolved,
        diffusive,
        total,
        turbulent,
        dispersive,
        True,
        sampled_resolved=sampled,
        snapshot_count=count,
    )


def read_turbine_data(
    config: CaseConfig,
    profile: ProfileData,
    start_time: float,
    end_time: float,
) -> TurbineData:
    paths = sorted(config.case_dir.glob("windfarm.*.nc"))
    if not paths:
        raise FileNotFoundError("No windfarm NetCDF output found")
    if len(paths) > 1:
        warnings.warn(f"Using first windfarm file: {paths[0].name}")

    rho_hub = float(np.interp(config.hub_height, profile.z, profile.rhoref))
    with nc.Dataset(paths[0]) as dataset:
        time = as_float_array(dataset.variables["time"])
        selected, _, _ = select_time_indices(time, start_time, end_time)
        weights = trapezoid_weights(time[selected])
        ids = np.asarray(dataset.variables["id"][:], dtype=int)
        filtered_velocity = as_float_array(dataset.variables["filtered_velocity"])
        thrust = as_float_array(dataset.variables["thrust"])

    plan_area = config.lx * config.ly / ids.size
    power_density = thrust * filtered_velocity / (rho_hub * plan_area)
    return TurbineData(
        time,
        ids,
        filtered_velocity,
        thrust,
        power_density,
        selected,
        weights,
        rho_hub,
        plan_area,
    )


def load_turbine_locations(config: CaseConfig) -> np.ndarray:
    locations = np.loadtxt(config.case_dir / "turbine_locations.txt")
    if locations.ndim != 2 or locations.shape[1] != 3:
        raise ValueError("Expected turbine location columns: id, x, y")
    return locations


def set_plot_style() -> None:
    plt.rcParams.update(
        {
            "font.size": 10,
            "axes.labelsize": 11,
            "axes.titlesize": 11,
            "legend.fontsize": 8.5,
            "figure.dpi": 120,
            "savefig.bbox": "tight",
            "axes.grid": True,
            "grid.alpha": 0.22,
            "grid.linewidth": 0.6,
            "lines.linewidth": 1.7,
        }
    )


def rotor_bounds(ax: plt.Axes, config: CaseConfig) -> None:
    ax.axhline(config.turbine_bottom / config.height, color="0.35", ls=":", lw=1.0)
    ax.axhline(config.turbine_top / config.height, color="0.35", ls=":", lw=1.0)


def save_figure(fig: plt.Figure, output: Path, dpi: int) -> None:
    fig.savefig(output, dpi=dpi)
    if not SHOW_FIGURES:
        plt.close(fig)
    print(f"Wrote {output.name}")


def nearest_time_index(time: np.ndarray, requested: float) -> int:
    return int(np.argmin(np.abs(time - requested)))


def plot_instantaneous_velocity(
    config: CaseConfig,
    locations: np.ndarray,
    end_time: float,
    output_dir: Path,
    dpi: int,
) -> None:
    paths = {
        "xz": config.case_dir / "u.xz.nc",
        "yz": config.case_dir / "u.yz.nc",
        "xy": config.case_dir / "u.xy.nc",
    }
    if any(not path.exists() for path in paths.values()):
        warnings.warn("Skipping Figure 1: one or more u cross-section files are missing")
        return

    panels: dict[str, tuple[np.ndarray, np.ndarray, np.ndarray, float]] = {}
    with nc.Dataset(paths["xz"]) as ds:
        time = as_float_array(ds.variables["time"])
        record = nearest_time_index(time, end_time)
        panels["xz"] = (
            as_float_array(ds.variables["xh"]),
            as_float_array(ds.variables["z"]),
            np.asarray(ds.variables["u"][record, :, 0, :], dtype=float),
            float(time[record]),
        )
    with nc.Dataset(paths["yz"]) as ds:
        time = as_float_array(ds.variables["time"])
        record = nearest_time_index(time, end_time)
        panels["yz"] = (
            as_float_array(ds.variables["y"]),
            as_float_array(ds.variables["z"]),
            np.asarray(ds.variables["u"][record, :, :, 0], dtype=float),
            float(time[record]),
        )
    with nc.Dataset(paths["xy"]) as ds:
        time = as_float_array(ds.variables["time"])
        record = nearest_time_index(time, end_time)
        panels["xy"] = (
            as_float_array(ds.variables["xh"]),
            as_float_array(ds.variables["y"]),
            np.asarray(ds.variables["u"][record, 0, :, :], dtype=float),
            float(time[record]),
        )

    fig, axes = plt.subplots(3, 1, figsize=(8.0, 9.5), constrained_layout=True)
    labels = {
        "xz": (r"$x/H$", r"$z/H$", "(a) streamwise-vertical plane"),
        "yz": (r"$y/H$", r"$z/H$", "(b) cross-stream vertical plane"),
        "xy": (r"$x/H$", r"$y/H$", "(c) hub-height horizontal plane"),
    }
    mesh = None
    for ax, mode in zip(axes, ("xz", "yz", "xy")):
        horizontal, vertical, velocity, time_value = panels[mode]
        ax.grid(False)
        mesh = ax.pcolormesh(
            horizontal / config.height,
            vertical / config.height,
            velocity / config.ustar,
            shading="auto",
            cmap="turbo",
            vmin=0.0,
            vmax=16.0,
            rasterized=True,
        )
        ax.set_xlabel(labels[mode][0])
        ax.set_ylabel(labels[mode][1])
        ax.set_title(f"{labels[mode][2]},  $t u_*/H={time_value*config.ustar/config.height:.2f}$")

        if mode == "xz":
            for x in np.unique(np.round(locations[:, 1], 3)):
                ax.plot(
                    [x / config.height, x / config.height],
                    [config.turbine_bottom / config.height, config.turbine_top / config.height],
                    color="k",
                    lw=2.0,
                )
        elif mode == "yz":
            for y in np.unique(np.round(locations[:, 2], 2)):
                ax.plot(
                    [y / config.height, y / config.height],
                    [config.turbine_bottom / config.height, config.turbine_top / config.height],
                    color="k",
                    lw=2.0,
                )
        else:
            for _, x, y in locations:
                ax.plot(
                    [x / config.height, x / config.height],
                    [(y - config.diameter / 2) / config.height, (y + config.diameter / 2) / config.height],
                    color="k",
                    lw=1.5,
                )

    assert mesh is not None
    colorbar = fig.colorbar(mesh, ax=axes, pad=0.02, aspect=35)
    colorbar.set_label(r"$u/u_*$")
    save_figure(fig, output_dir / "figure01_instantaneous_velocity.png", dpi)


def plot_mean_velocity(
    config: CaseConfig,
    profile: ProfileData,
    means: dict[str, np.ndarray | float],
    output_dir: Path,
    dpi: int,
) -> None:
    mean_u = np.asarray(means["u"])
    z_normalized = profile.z / config.hub_height
    log_reference = np.log(profile.z / config.roughness) / KAPPA

    fig, ax = plt.subplots(figsize=(6.0, 5.3), constrained_layout=True)
    ax.semilogx(z_normalized, mean_u / config.ustar, color="k", label="MicroHH")
    ax.semilogx(
        z_normalized,
        log_reference,
        color="0.55",
        ls="--",
        label=r"no-turbine log law",
    )
    ax.axvline(config.turbine_bottom / config.hub_height, color="0.35", ls=":")
    ax.axvline(config.turbine_top / config.hub_height, color="0.35", ls=":")
    ax.set_xlabel(r"$z/z_h$")
    ax.set_ylabel(r"$\langle\overline{u}\rangle/u_*$")
    ax.set_title("Mean velocity profile (Calaf Figs. 2 and 7)")
    ax.set_xlim(max(0.03, z_normalized[0]), z_normalized[-1])
    ax.set_ylim(bottom=0.0)
    ax.legend()
    save_figure(fig, output_dir / "figure02_mean_velocity.png", dpi)


def mean_velocity_at_half_levels(profile: ProfileData, mean_u: np.ndarray) -> np.ndarray:
    return np.interp(profile.zh, profile.z, mean_u, left=mean_u[0], right=mean_u[-1])


def mean_gradient_at_half_levels(profile: ProfileData, mean_u: np.ndarray) -> np.ndarray:
    gradient_full = np.gradient(mean_u, profile.z, edge_order=2)
    return np.interp(
        profile.zh,
        profile.z,
        gradient_full,
        left=gradient_full[0],
        right=gradient_full[-1],
    )


def plot_shear_stress(
    config: CaseConfig,
    profile: ProfileData,
    stresses: StressData,
    output_dir: Path,
    dpi: int,
) -> None:
    scale = config.ustar**2
    y = profile.zh / config.height
    fig, ax = plt.subplots(figsize=(6.2, 5.3), constrained_layout=True)

    if stresses.separation_available:
        ax.plot(stresses.turbulent / scale, y, ls="-.", label="turbulent resolved")
        ax.plot(stresses.dispersive / scale, y, ls="--", label="dispersive")
        ax.plot(stresses.resolved / scale, y, color="k", label="resolved sum")
    else:
        ax.plot(stresses.resolved / scale, y, color="k", label="resolved (not separated)")
    ax.plot(stresses.diffusive / scale, y, ls=":", label="SGS + molecular")
    ax.plot(stresses.total / scale, y, color="tab:red", lw=1.2, label="total incl. SGS")
    rotor_bounds(ax, config)
    ax.set_xlabel(r"$\tau_{xz}/u_*^2$")
    ax.set_ylabel(r"$z/H$")
    ax.set_title("Vertical shear-stress profiles (Calaf Figs. 3 and 8)")
    ax.set_ylim(0.0, 1.0)
    ax.legend(loc="best")
    save_figure(fig, output_dir / "figure03_shear_stress.png", dpi)


def plot_energy_conversion(
    config: CaseConfig,
    profile: ProfileData,
    mean_u: np.ndarray,
    stresses: StressData,
    output_dir: Path,
    dpi: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    gradient = mean_gradient_at_half_levels(profile, mean_u)
    normalization = config.ustar**3 / config.height
    turbulent = stresses.turbulent * gradient
    dispersive = stresses.dispersive * gradient
    resolved = stresses.resolved * gradient

    fig, ax = plt.subplots(figsize=(6.2, 5.3), constrained_layout=True)
    if stresses.separation_available:
        ax.plot(turbulent / normalization, profile.zh / config.height, ls="-.", label="turbulent")
        ax.plot(dispersive / normalization, profile.zh / config.height, ls="--", label="dispersive")
    ax.plot(resolved / normalization, profile.zh / config.height, color="k", label="resolved sum")
    rotor_bounds(ax, config)
    ax.axvline(0.0, color="0.55", lw=0.8)
    ax.set_xlabel(r"$\tau_{xz}\,\partial_z U/(u_*^3/H)$")
    ax.set_ylabel(r"$z/H$")
    ax.set_title("Mean-energy conversion (Calaf Fig. 4)")
    ax.set_ylim(0.0, 1.0)
    ax.legend()
    save_figure(fig, output_dir / "figure04_energy_conversion.png", dpi)
    return turbulent, dispersive, resolved


def plot_energy_flux(
    config: CaseConfig,
    profile: ProfileData,
    mean_u: np.ndarray,
    stresses: StressData,
    output_dir: Path,
    dpi: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    u_half = mean_velocity_at_half_levels(profile, mean_u)
    turbulent = stresses.turbulent * u_half
    dispersive = stresses.dispersive * u_half
    resolved = stresses.resolved * u_half

    fig, ax = plt.subplots(figsize=(6.2, 5.3), constrained_layout=True)
    if stresses.separation_available:
        ax.plot(turbulent / config.ustar**3, profile.zh / config.height, ls="-.", label="turbulent")
        ax.plot(dispersive / config.ustar**3, profile.zh / config.height, ls="--", label="dispersive")
    ax.plot(resolved / config.ustar**3, profile.zh / config.height, color="k", label="resolved sum")
    rotor_bounds(ax, config)
    ax.axvline(0.0, color="0.55", lw=0.8)
    ax.set_xlabel(r"$\Phi/u_*^3$")
    ax.set_ylabel(r"$z/H$")
    ax.set_title("Vertical kinetic-energy flux (Calaf Fig. 5)")
    ax.set_ylim(0.0, 1.0)
    ax.legend()
    save_figure(fig, output_dir / "figure05_energy_flux.png", dpi)
    return turbulent, dispersive, resolved


def group_turbines(locations: np.ndarray) -> tuple[list[int], list[int]]:
    # Round away sub-millimetre differences in the hand-written location file.
    rounded_x = np.round(locations[:, 1], 2)
    rounded_y = np.round(locations[:, 2], 2)
    first_y = np.min(rounded_y)
    first_x = np.min(rounded_x)
    same_row = locations[np.isclose(rounded_y, first_y)]
    same_column = locations[np.isclose(rounded_x, first_x)]
    same_row = same_row[np.argsort(same_row[:, 1])]
    same_column = same_column[np.argsort(same_column[:, 2])]
    return same_row[:, 0].astype(int).tolist(), same_column[:4, 0].astype(int).tolist()


def plot_turbine_power(
    config: CaseConfig,
    turbine: TurbineData,
    locations: np.ndarray,
    start_time: float,
    end_time: float,
    output_dir: Path,
    dpi: int,
) -> None:
    same_row, same_column = group_turbines(locations)
    id_to_index = {int(turbine_id): i for i, turbine_id in enumerate(turbine.ids)}
    nondimensional_time = turbine.time * config.ustar / config.height

    fig, axes = plt.subplots(2, 1, figsize=(8.0, 6.5), sharex=True, constrained_layout=True)
    for ax, ids, title in (
        (axes[0], same_row, "(a) same row, different streamwise columns"),
        (axes[1], same_column, "(b) same streamwise column, different rows"),
    ):
        for turbine_id in ids:
            index = id_to_index[turbine_id]
            ax.plot(
                nondimensional_time,
                turbine.power_density[:, index] / config.ustar**3,
                label=f"ID {turbine_id}",
                lw=1.0,
            )
        ax.axvspan(
            start_time * config.ustar / config.height,
            end_time * config.ustar / config.height,
            color="0.85",
            alpha=0.35,
            zorder=-10,
            label="averaging window" if ax is axes[0] else None,
        )
        ax.set_ylabel(r"$P_k/u_*^3$")
        ax.set_title(title)
        ax.legend(ncol=3, loc="best")
    axes[-1].set_xlabel(r"$t u_*/H$")
    fig.suptitle("Extracted turbine power density (Calaf Fig. 6)")
    save_figure(fig, output_dir / "figure06_turbine_power.png", dpi)


def plot_convergence(
    config: CaseConfig,
    profile: ProfileData,
    turbine: TurbineData,
    start_time: float,
    end_time: float,
    output_dir: Path,
    dpi: int,
) -> None:
    u_hub = np.array(
        [np.interp(config.hub_height, profile.z, values) for values in profile.u]
    )
    turbine_mean = np.mean(turbine.power_density, axis=1)
    fig, axes = plt.subplots(3, 1, figsize=(8.0, 7.2), sharex=True, constrained_layout=True)
    axes[0].plot(profile.time * config.ustar / config.height, u_hub / config.ustar)
    axes[0].set_ylabel(r"$U_h/u_*$")
    axes[1].plot(
        profile.time * config.ustar / config.height,
        profile.ustar_surface / config.ustar,
    )
    axes[1].set_ylabel(r"$u_{*,lo}/u_*$")
    axes[2].plot(turbine.time * config.ustar / config.height, turbine_mean / config.ustar**3)
    axes[2].set_ylabel(r"$\langle P_k\rangle/u_*^3$")
    axes[2].set_xlabel(r"$t u_*/H$")
    for ax in axes:
        ax.axvspan(
            start_time * config.ustar / config.height,
            end_time * config.ustar / config.height,
            color="tab:blue",
            alpha=0.10,
        )
    fig.suptitle("Convergence indicators (shading marks averaging window)")
    save_figure(fig, output_dir / "convergence.png", dpi)


def interpolate_finite(x_new: float, x: np.ndarray, values: np.ndarray) -> float:
    valid = np.isfinite(x) & np.isfinite(values)
    if np.count_nonzero(valid) < 2:
        return float("nan")
    return float(np.interp(x_new, x[valid], values[valid]))


def layer_average(z: np.ndarray, values: np.ndarray, bottom: float, top: float) -> float:
    valid = np.isfinite(z) & np.isfinite(values)
    z_valid = z[valid]
    values_valid = values[valid]
    interior = (z_valid > bottom) & (z_valid < top)
    sample_z = np.concatenate(([bottom], z_valid[interior], [top]))
    sample_values = np.interp(sample_z, z_valid, values_valid)
    return float(np.trapz(sample_values, sample_z) / (top - bottom))


def layer_integral(z: np.ndarray, values: np.ndarray, bottom: float, top: float) -> float:
    valid = np.isfinite(z) & np.isfinite(values)
    z_valid = z[valid]
    values_valid = values[valid]
    interior = (z_valid > bottom) & (z_valid < top)
    sample_z = np.concatenate(([bottom], z_valid[interior], [top]))
    sample_values = np.interp(sample_z, z_valid, values_valid)
    return float(np.trapz(sample_values, sample_z))


def write_csv_outputs(
    config: CaseConfig,
    profile: ProfileData,
    means: dict[str, np.ndarray | float],
    stresses: StressData,
    production: tuple[np.ndarray, np.ndarray, np.ndarray],
    energy_flux: tuple[np.ndarray, np.ndarray, np.ndarray],
    output_dir: Path,
) -> None:
    mean_u = np.asarray(means["u"])
    mean_table = np.column_stack(
        (
            profile.z,
            profile.z / config.hub_height,
            mean_u,
            mean_u / config.ustar,
            np.sqrt(np.maximum(np.asarray(means["u2"]), 0.0)),
        )
    )
    np.savetxt(
        output_dir / "mean_velocity_profile.csv",
        mean_table,
        delimiter=",",
        header="z_m,z_over_zh,U_m_s,U_over_ustar,instantaneous_plane_rms_u_m_s",
        comments="",
    )

    prod_turb, prod_disp, prod_res = production
    flux_turb, flux_disp, flux_res = energy_flux
    scale2 = config.ustar**2
    scale3 = config.ustar**3
    flux_table = np.column_stack(
        (
            profile.zh,
            profile.zh / config.height,
            stresses.turbulent / scale2,
            stresses.dispersive / scale2,
            stresses.resolved / scale2,
            stresses.diffusive / scale2,
            stresses.total / scale2,
            prod_turb / (scale3 / config.height),
            prod_disp / (scale3 / config.height),
            prod_res / (scale3 / config.height),
            flux_turb / scale3,
            flux_disp / scale3,
            flux_res / scale3,
        )
    )
    np.savetxt(
        output_dir / "stress_energy_profiles.csv",
        flux_table,
        delimiter=",",
        header=(
            "zh_m,z_over_H,tau_turb_over_ustar2,tau_disp_over_ustar2,"
            "tau_resolved_over_ustar2,tau_sgs_over_ustar2,tau_total_over_ustar2,"
            "production_turb_normalized,production_disp_normalized,"
            "production_resolved_normalized,flux_turb_over_ustar3,"
            "flux_disp_over_ustar3,flux_resolved_over_ustar3"
        ),
        comments="",
    )
    print("Wrote mean_velocity_profile.csv")
    print("Wrote stress_energy_profiles.csv")


def write_summary(
    config: CaseConfig,
    profile: ProfileData,
    means: dict[str, np.ndarray | float],
    stresses: StressData,
    turbine: TurbineData,
    production: tuple[np.ndarray, np.ndarray, np.ndarray],
    energy_flux: tuple[np.ndarray, np.ndarray, np.ndarray],
    start_time: float,
    end_time: float,
    output_dir: Path,
) -> None:
    del production  # The resolved integral is recomputed below for clarity.
    mean_u = np.asarray(means["u"])
    u_hub = float(np.interp(config.hub_height, profile.z, mean_u))
    u_2hub = float(np.interp(2.0 * config.hub_height, profile.z, mean_u))
    z0_hi = 2.0 * config.hub_height * np.exp(-KAPPA * u_2hub / config.ustar)
    u_disk_layer = layer_average(
        profile.z, mean_u, config.turbine_bottom, config.turbine_top
    )

    tau_res_bottom = interpolate_finite(
        config.turbine_bottom, profile.zh, stresses.resolved
    )
    tau_total_bottom = interpolate_finite(
        config.turbine_bottom, profile.zh, stresses.total
    )
    ustar_lo_resolved = np.sqrt(max(tau_res_bottom, 0.0))
    ustar_lo_total = np.sqrt(max(tau_total_bottom, 0.0))

    power_selected = turbine.power_density[turbine.selected]
    mean_power_by_turbine = weighted_mean(power_selected, turbine.weights)
    pt = float(np.mean(mean_power_by_turbine))
    mean_thrust_by_turbine = weighted_mean(
        turbine.thrust[turbine.selected], turbine.weights
    )
    mean_thrust = float(np.mean(mean_thrust_by_turbine))

    _, _, resolved_energy_flux = energy_flux
    phi_bottom = interpolate_finite(
        config.turbine_bottom, profile.zh, resolved_energy_flux
    )
    phi_top = interpolate_finite(config.turbine_top, profile.zh, resolved_energy_flux)
    delta_phi = phi_top - phi_bottom
    dudz = mean_gradient_at_half_levels(profile, mean_u)
    dissipation = layer_integral(
        profile.zh,
        stresses.resolved * dudz,
        config.turbine_bottom,
        config.turbine_top,
    )
    pressure_work = config.diameter * (-config.dpdx) * u_disk_layer
    turbine_mean_work = mean_thrust * u_disk_layer / (
        turbine.rho_hub * turbine.plan_area_per_turbine
    )
    scale3 = config.ustar**3
    budget_residual = (
        pressure_work + delta_phi - dissipation - turbine_mean_work
    ) / scale3

    duration_star = (end_time - start_time) * config.ustar / config.height
    total_run_star = profile.time[-1] * config.ustar / config.height
    lines = [
        "Calaf et al. (2010) analysis of MicroHH calaf2010",
        "=" * 58,
        "",
        f"Averaging window: {start_time:.3f} to {end_time:.3f} s",
        f"Averaging duration: {duration_star:.4f} H/u_*",
        f"Total simulated duration: {total_run_star:.4f} H/u_*",
        f"Statistics records: {profile.selected.size}",
        f"3-D snapshots used: {stresses.snapshot_count}",
        "",
        "Velocity and roughness",
        "----------------------",
        f"Pressure-gradient u_*: {config.ustar:.6f} m s-1",
        f"Mean surface u_*lo: {float(means['ustar_surface']):.6f} m s-1",
        f"Surface u_*lo/u_*: {float(means['ustar_surface'])/config.ustar:.6f}",
        f"Hub-height U/u_*: {u_hub/config.ustar:.6f}",
        f"Rotor-layer U_D/u_*: {u_disk_layer/config.ustar:.6f}",
        f"Effective z0_hi: {z0_hi:.6f} m",
        f"Effective z0_hi/z_h: {z0_hi/config.hub_height:.6f}",
        f"Stress-derived u_*lo^tu/u_* (resolved): {ustar_lo_resolved/config.ustar:.6f}",
        f"Stress-derived u_*lo/u_* (incl. SGS): {ustar_lo_total/config.ustar:.6f}",
        "Paper A1 reference: z0_hi/z_h=0.034, wall u_*lo/u_*=0.58, "
        "u_*lo^tu/u_*=0.55",
        "",
        "Rotor-layer mean-energy budget (normalized by u_*^3)",
        "----------------------------------------------------",
        f"Pressure work W_p/u_*^3: {pressure_work/scale3:.6f}",
        f"Vertical flux difference deltaPhi/u_*^3: {delta_phi/scale3:.6f}",
        f"Resolved conversion D/u_*^3: {dissipation/scale3:.6f}",
        f"Mean-flow turbine work W_T/u_*^3: {turbine_mean_work/scale3:.6f}",
        f"Budget residual (W_p+deltaPhi-D-W_T)/u_*^3: {budget_residual:.6f}",
        f"Disk-velocity turbine power P_T/u_*^3: {pt/scale3:.6f}",
        "Paper A2 references: W_p=0.85, deltaPhi=5.4, D=0.68, W_T=5.7, P_T=4.93",
        "",
        "Interpretation notes",
        "--------------------",
        "- Calaf used 60 H/u_* of startup before the A1 averaging period.",
        "- The present short run is not expected to be statistically stationary.",
        "- Dispersive stress is estimated from the 10-minute 3-D snapshots.",
        "- Turbulent stress is resolved total stress minus dispersive stress.",
        "- SGS stress is reported separately and is not included in Calaf's resolved D and Phi terms.",
        "- MicroHH ADM writes power=0; P_T here is thrust times filtered disk velocity.",
    ]

    if stresses.sampled_resolved is not None:
        valid = (
            np.isfinite(stresses.sampled_resolved)
            & np.isfinite(stresses.resolved)
            & (profile.zh >= config.turbine_bottom)
            & (profile.zh <= config.turbine_top)
        )
        if np.any(valid):
            rms = np.sqrt(
                np.mean(
                    (stresses.sampled_resolved[valid] - stresses.resolved[valid]) ** 2
                )
            )
            lines.append(
                f"- 3-D-snapshot versus 1-minute resolved-stress RMS difference: {rms/config.ustar**2:.6f} u_*^2."
            )

    summary_path = output_dir / "calaf2010_summary.txt"
    summary_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"Wrote {summary_path.name}")


def main() -> None:
    global SHOW_FIGURES

    args = parse_arguments()
    SHOW_FIGURES = args.show
    case_dir = args.case_dir.resolve()
    config = read_case_config(case_dir)
    output_dir = (
        args.output_dir.resolve()
        if args.output_dir is not None
        else case_dir / "calaf2010_figures"
    )
    output_dir.mkdir(parents=True, exist_ok=True)
    set_plot_style()

    profile, start_time, end_time = read_profile_data(
        config, args.start_time, args.end_time
    )
    means = time_mean_profiles(profile)
    duration_star = (end_time - start_time) * config.ustar / config.height
    print(
        f"Averaging {start_time:.1f} <= t <= {end_time:.1f} s "
        f"({duration_star:.3f} H/u_*)"
    )
    if duration_star < 60.0:
        warnings.warn(
            f"Averaging interval is only {duration_star:.2f} H/u_*; "
            "Calaf A1 accumulated statistics for 60 H/u_* after spinup."
        )

    stresses = assemble_stresses(
        profile, means, config, start_time, end_time, args.skip_3d
    )
    turbine = read_turbine_data(config, profile, start_time, end_time)
    locations = load_turbine_locations(config)

    if not args.skip_instantaneous:
        plot_instantaneous_velocity(
            config, locations, end_time, output_dir, args.dpi
        )
    plot_mean_velocity(config, profile, means, output_dir, args.dpi)
    plot_shear_stress(config, profile, stresses, output_dir, args.dpi)
    production = plot_energy_conversion(
        config, profile, np.asarray(means["u"]), stresses, output_dir, args.dpi
    )
    energy_flux = plot_energy_flux(
        config, profile, np.asarray(means["u"]), stresses, output_dir, args.dpi
    )
    plot_turbine_power(
        config,
        turbine,
        locations,
        start_time,
        end_time,
        output_dir,
        args.dpi,
    )
    plot_convergence(
        config,
        profile,
        turbine,
        start_time,
        end_time,
        output_dir,
        args.dpi,
    )
    write_csv_outputs(
        config, profile, means, stresses, production, energy_flux, output_dir
    )
    write_summary(
        config,
        profile,
        means,
        stresses,
        turbine,
        production,
        energy_flux,
        start_time,
        end_time,
        output_dir,
    )

    print(f"Analysis complete: {output_dir}")
    if args.show:
        plt.show()


if __name__ == "__main__":
    main()
