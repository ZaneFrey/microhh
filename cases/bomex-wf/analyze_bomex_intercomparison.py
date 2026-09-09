#!/usr/bin/env python3
"""Create BOMEX intercomparison-style diagnostics from MicroHH statistics.

The script uses the 3--6 h averaging interval of Siebesma et al. (2003).
It expects the three statistics files written by a BOMEX run in the case
directory and writes PNG files to ``outputs/``.  Figures 9, 10, and the
right-side diagnostic in Figure 12 are LES-derived approximations: the
published curves additionally use a bulk plume model and cloud schemes that
are not part of the MicroHH statistics output.

Run from this directory (or provide --case-dir):
    python3 analyze_bomex_intercomparison.py
"""

from __future__ import annotations

import argparse
from pathlib import Path
import warnings

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from netCDF4 import Dataset


AVG_START = 3 * 3600.0
AVG_END = 6 * 3600.0
KM = 1.0e-3
G = 9.81
REF_Z = np.array([0.0, 0.5, 1.0, 1.5, 2.0, 2.5])


def style() -> None:
    plt.rcParams.update({
        "figure.dpi": 140, "savefig.dpi": 180, "font.size": 10,
        "axes.grid": True, "grid.alpha": 0.25, "axes.spines.top": False,
        "axes.spines.right": False, "lines.linewidth": 2.0,
    })


def load(path: Path) -> dict:
    """Read a grouped MicroHH statistics NetCDF file into a simple dict."""
    with Dataset(path) as nc:
        out = {name: np.asarray(nc[name][:]) for name in ("time", "z", "zh") if name in nc.variables}
        for group_name, group in nc.groups.items():
            for name, var in group.variables.items():
                values = var[:]
                out[f"{group_name}/{name}"] = np.asarray(
                    values.filled(np.nan) if np.ma.isMaskedArray(values) else values
                )
    return out


def load_restart_series(case: Path, basename: str) -> dict:
    """Load and join statistics written across one or more restarts.

    Restarted MicroHH runs write a new NetCDF file whose suffix is the
    restart time.  Failed attempts can overlap both the original file and a
    later restart.  Files are applied in increasing restart-time order and
    the newest segment wins wherever sample times overlap.
    """
    paths = sorted(
        case.glob(f"{basename}.[0-9][0-9][0-9][0-9][0-9][0-9][0-9].nc"),
        key=lambda path: int(path.stem.rsplit(".", 1)[1]),
    )
    if not paths:
        raise FileNotFoundError(f"No statistics files found for {basename}")

    base_paths = [path for path in paths if path.stem.endswith(".0000000")]
    if base_paths:
        base_mtime = base_paths[0].stat().st_mtime
        stale_paths = [path for path in paths if path.stat().st_mtime < base_mtime]
        if stale_paths:
            warnings.warn(
                "Ignoring restart statistics older than the current base file: "
                + ", ".join(path.name for path in stale_paths)
            )
            paths = [path for path in paths if path not in stale_paths]

    segments = [load(path) for path in paths]
    for coordinate in ("z", "zh"):
        reference = segments[0][coordinate]
        if any(
            segment[coordinate].shape != reference.shape
            or not np.allclose(segment[coordinate], reference)
            for segment in segments[1:]
        ):
            raise ValueError(f"Coordinate {coordinate} changes between {basename} restart files")

    all_times = np.concatenate([segment["time"] for segment in segments])
    # Adaptive stepping can leave nominal output times differing by roughly
    # 1e-10 s between restart attempts.  Quantize to microseconds before
    # identifying overlaps so 14460.0 and 14459.999999999995 are one sample.
    time_keys = np.rint(all_times * 1.0e6).astype(np.int64)
    # np.unique keeps its first match.  Reverse first so that duplicate times
    # are taken from the later restart segment, then restore chronological order.
    _, reverse_indices = np.unique(time_keys[::-1], return_index=True)
    keep = all_times.size - 1 - reverse_indices
    keep = keep[np.argsort(all_times[keep])]

    merged = {
        "time": all_times[keep],
        "z": segments[-1]["z"],
        "zh": segments[-1]["zh"],
    }
    keys = set.intersection(*(set(segment) for segment in segments)) - {"time", "z", "zh"}
    for key in keys:
        values = [segment[key] for segment in segments]
        if all(value.ndim > 0 and value.shape[0] == segment["time"].size
               for value, segment in zip(values, segments)):
            merged[key] = np.concatenate(values, axis=0)[keep]
        else:
            merged[key] = values[-1]

    source_names = ", ".join(path.name for path in paths)
    print(f"Loaded {basename} from {source_names} ({merged['time'].size} unique samples)")
    return merged


def field(data: dict, group: str, name: str) -> np.ndarray:
    key = f"{group}/{name}"
    if key not in data:
        raise KeyError(f"Required statistic {key} is missing")
    return data[key]


def mean_last3(data: dict, group: str, name: str) -> np.ndarray:
    times = data["time"]
    use = (times >= AVG_START) & (times <= AVG_END)
    if not np.any(use):
        raise ValueError("No samples in the 3--6 h averaging interval")
    # Conditional fields are intentionally undefined where a mask has no
    # points (for example above cloud top).  Preserve NaN there without a
    # distracting NumPy warning.
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        return np.nanmean(field(data, group, name)[use], axis=0)


def save(fig: plt.Figure, output: Path, name: str) -> None:
    fig.tight_layout()
    fig.savefig(output / name, bbox_inches="tight")
    plt.close(fig)


def ref_band(ax: plt.Axes, lo: np.ndarray, hi: np.ndarray, z: np.ndarray = REF_Z,
             label: str | None = None) -> None:
    """Plot a digitized Siebesma et al. ensemble envelope behind LES data.

    Values are graphical read-offs of the paper's gray, +/-2-standard-
    deviation bands at the tabulated heights; interpolation gives a smooth
    envelope.  They are deliberately kept separate from this run's sampling
    uncertainty.
    """
    ax.fill_betweenx(z, lo, hi, color="0.5", alpha=0.50, lw=0, zorder=0,
                     label=label)


def digitized_band(ax: plt.Axes, reference: np.lib.npyio.NpzFile, name: str, scale: float = 1.0) -> None:
    """Draw the pixel-level envelope digitized from a supplied paper crop."""
    z = reference[f"{name}_z"] * KM
    ax.fill_betweenx(z, reference[f"{name}_lo"] * scale, reference[f"{name}_hi"] * scale,
                     color="0.5", alpha=0.50, lw=0, zorder=0)


def digitized_time_band(ax: plt.Axes, reference: np.lib.npyio.NpzFile, name: str, scale: float = 1.0) -> None:
    ax.fill_between(reference[f"{name}_x"] / 60., reference[f"{name}_lo"] * scale,
                    reference[f"{name}_hi"] * scale, color="0.5", alpha=0.50, lw=0, zorder=0)


def prof_axis(ax: plt.Axes) -> None:
    ax.set_ylim(0, 3.5)
    ax.set_ylabel("height (km)")


def safe_div(num: np.ndarray, den: np.ndarray) -> np.ndarray:
    result = np.full_like(np.asarray(num, dtype=float), np.nan)
    good = np.abs(den) > 1.0e-14
    result[good] = np.asarray(num, dtype=float)[good] / np.asarray(den, dtype=float)[good]
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case-dir", type=Path, default=Path(__file__).resolve().parent)
    args = parser.parse_args()
    case = args.case_dir.resolve()
    output = case / "outputs"
    output.mkdir(exist_ok=True)

    default = load_restart_series(case, "bomex-wf.default")
    cloud = load_restart_series(case, "bomex-wf.ql")
    core = load_restart_series(case, "bomex-wf.qlcore")
    initial = load(case / "bomex-wf_input.nc")
    reference_path = case / "siebesma_envelopes.npz"
    if not reference_path.exists():
        archived_reference = case / "old" / "spinup" / reference_path.name
        if archived_reference.exists():
            reference_path = archived_reference
        else:
            raise FileNotFoundError("Run digitize_siebesma_envelopes.py before this analysis.")
    reference = np.load(reference_path)
    z, zh = default["z"], default["zh"]
    zk, zhk = z * KM, zh * KM
    wh_to_z_2d = lambda values: 0.5 * (values[:, :-1] + values[:, 1:])
    style()

    # Figure 1 -- initial state.
    fig, axs = plt.subplots(1, 3, figsize=(10, 4.5), sharey=True)
    axs[0].plot(field(initial, "init", "thl"), zk, color="black"); axs[0].set(xlabel=r"$\theta_l$ (K)")
    axs[1].plot(1e3 * field(initial, "init", "qt"), zk, color="black"); axs[1].set(xlabel=r"$q_t$ (g kg$^{-1}$)")
    axs[2].plot(field(initial, "init", "u"), zk, color="black", label="u")
    axs[2].plot(field(initial, "init", "v"), zk, color="black", ls="--", label="v")
    axs[2].set(xlabel="wind (m s$^{-1}$)"); axs[2].legend()
    for ax in axs: prof_axis(ax)
    fig.suptitle("BOMEX Figure 1 - initial profiles")
    save(fig, output, "bomex_figure_01_initial_profiles.png")

    # Figure 2 -- time evolution.
    time_h = default["time"] / 3600.0
    # ``tke`` is present in the file, but reconstruct from variances so this
    # remains valid for MicroHH versions that leave its diagnostic unfilled.
    tke_all = 0.5 * (field(default, "default", "u_2") + field(default, "default", "v_2")
                     + wh_to_z_2d(field(default, "default", "w_2")))
    tke_column = np.trapezoid(tke_all, z, axis=1)
    fig, axs = plt.subplots(1, 3, figsize=(11, 3.6), sharex=True)
    digitized_time_band(axs[0], reference, "f2_cover", 100)
    digitized_time_band(axs[1], reference, "f2_lwp")
    digitized_time_band(axs[2], reference, "f2_tke")
    axs[0].plot(time_h, 100 * field(default, "thermo", "ql_cover"), color="black"); axs[0].set(ylabel="cloud cover (%)")
    axs[1].plot(time_h, 1e3 * field(default, "thermo", "ql_path"), color="black"); axs[1].set(ylabel="LWP (g m$^{-2}$)")
    axs[2].plot(time_h, tke_column, color="black"); axs[2].set(ylabel=r"$\int$ TKE dz (m$^3$ s$^{-2}$)")
    for ax in axs:
        ax.axvspan(3, 6, color="0.85", zorder=0, label="analysis interval")
        ax.set(xlim=(0, 6), xlabel="time (h)")
    axs[0].legend(loc="best")
    fig.suptitle("BOMEX Figure 2 - domain evolution")
    save(fig, output, "bomex_figure_02_time_series.png")

    # Convenience profiles from the final three hours.
    d = lambda group, name: mean_last3(default, group, name)
    c = lambda group, name: mean_last3(cloud, group, name)
    co = lambda group, name: mean_last3(core, group, name)
    wh_to_z = lambda values: 0.5 * (values[:-1] + values[1:])
    thl, qt, ql, thv = (d("thermo", x) for x in ("thl", "qt", "ql", "thv"))
    u, v = d("default", "u"), d("default", "v")

    # Figure 3 -- mean profiles.
    fig, axs = plt.subplots(1, 5, figsize=(14, 4.5), sharey=True)
    panels = [(thl, field(initial, "init", "thl"), r"$\theta_l$ (K)", "f3_thl"),
              (1e3 * qt, 1e3 * field(initial, "init", "qt"), r"$q_t$ (g kg$^{-1}$)", "f3_qt"),
              (u, field(initial, "init", "u"), "u (m s$^{-1}$)", "f3_u"),
              (v, field(initial, "init", "v"), "v (m s$^{-1}$)", None),
              (1e3 * ql, np.zeros_like(z), r"$q_l$ (g kg$^{-1}$)", "f3_ql")]
    for ax, (values, initial_values, label, reference_name) in zip(axs, panels):
        if reference_name:
            digitized_band(ax, reference, reference_name)
        ax.plot(values, zk, color="black")
        ax.plot(initial_values, zk, color="black", ls="--", lw=1.4)
        ax.set_xlabel(label); prof_axis(ax)
    fig.suptitle("BOMEX Figure 3 - mean profiles, 3-6 h")
    save(fig, output, "bomex_figure_03_mean_profiles.png")

    # Figure 4 -- total turbulent fluxes.
    fig, axs = plt.subplots(1, 5, figsize=(14, 4.5), sharey=True)
    # The paper labels scalar fluxes in W m-2.  These envelopes are converted
    # to native cloud-statistics units using rho=1.2 kg m-3, cp=1004 J kg-1
    # K-1, and Lv=2.5e6 J kg-1; no conversion is applied to this LES output.
    fluxes = [(1e3*d("thermo", "qt_flux"), r"$w'q_t'$ (g kg$^{-1}$ m s$^{-1}$)", "f4_qt", 1e3/(1.2*2.5e6)),
              (d("thermo", "thl_flux"), r"$w'\theta_l'$ (K m s$^{-1}$)", "f4_thl", 1/(1.2*1004.)),
              (1e3*d("thermo", "ql_flux"), r"$w'q_l'$ (g kg$^{-1}$ m s$^{-1}$)", "f4_ql", 1e3/(1.2*2.5e6)),
              (d("thermo", "thv_flux"), r"$w'\theta_v'$ (K m s$^{-1}$)", "f4_thv", 1/(1.2*1004.)),
              (d("default", "u_flux"), r"$u'w'$ (m$^2$ s$^{-2}$)", "f4_uw", 1)]
    for index, (ax, (values, label, reference_name, scale)) in enumerate(zip(axs, fluxes)):
        digitized_band(ax, reference, reference_name, scale)
        ax.plot(values, zhk, color="black"); ax.set_xlabel(label); prof_axis(ax)
        if index:
            ax.set_ylabel("")
    fig.suptitle("BOMEX Figure 4 - turbulent-flux profiles, 3-6 h")
    save(fig, output, "bomex_figure_04_flux_profiles.png")

    # Figure 5 -- turbulent kinetic energy and vertical velocity variance.
    fig, axs = plt.subplots(1, 2, figsize=(7.4, 4.5), sharey=True)
    use = (default["time"] >= AVG_START) & (default["time"] <= AVG_END)
    digitized_band(axs[0], reference, "f5_tke")
    digitized_band(axs[1], reference, "f5_w2")
    axs[0].plot(np.nanmean(tke_all[use], axis=0), zk, color="black"); axs[0].set(xlabel="TKE (m$^2$ s$^{-2}$)")
    axs[1].plot(d("default", "w_2"), zhk, color="black"); axs[1].set(xlabel=r"$w'^2$ (m$^2$ s$^{-2}$)")
    for ax in axs: prof_axis(ax)
    fig.suptitle("BOMEX Figure 5 - turbulence profiles, 3-6 h")
    save(fig, output, "bomex_figure_05_tke_wvariance.png")

    # Figure 6 -- cloud and positively buoyant cloud-core cover.
    cloud_area, core_area = c("default", "area"), co("default", "area")
    fig, ax = plt.subplots(figsize=(5, 4.5))
    digitized_band(ax, reference, "f6_cover", 100)
    ax.plot(100 * cloud_area, zk, label="cloud (ql > 0)")
    ax.plot(100 * core_area, zk, label="cloud core (ql > 0, b' > 0)")
    ax.set(xlabel="fractional area (%)"); prof_axis(ax); ax.legend()
    fig.suptitle("BOMEX Figure 6 - cloud and core cover, 3-6 h")
    save(fig, output, "bomex_figure_06_cloud_core_cover.png")

    # Figure 7 -- conditional cloud/core values.
    fig, axs = plt.subplots(1, 5, figsize=(14, 4.5), sharey=True)
    conditional = [("thl", r"$\theta_l$ (K)", 1), ("qt", r"$q_t$ (g kg$^{-1}$)", 1e3),
                   ("thv", r"$\theta_v$ (K)", 1), ("ql", r"$q_l$ (g kg$^{-1}$)", 1e3),
                   ("w", "w (m s$^{-1}$)", 1)]
    # Approximate cloud/core envelopes digitized from Fig. 7 (cloud first,
    # core second); values beyond the reliably sampled cloud layer are omitted.
    bands7 = {
        "thl": (([298.9, 299.2, 299.8, 300.5, 301.4, np.nan], [.10, .15, .25, .30, .30, np.nan]),
                ([298.9, 299.2, 299.7, 300.1, 300.8, np.nan], [.10, .12, .20, .25, .25, np.nan])),
        "qt": (([17.1, 16.6, 15.3, 14.0, 13.3, np.nan], [.2, .3, .5, .6, .4, np.nan]),
               ([17.1, 16.7, 15.6, 14.6, 14.2, np.nan], [.2, .25, .4, .5, .35, np.nan])),
        "thv": (([302.0, 302.0, 303.0, 304.7, 307.0, np.nan], [.1, .15, .30, .40, .35, np.nan]),
                ([302.0, 302.0, 302.8, 304.2, 306.0, np.nan], [.1, .12, .25, .35, .35, np.nan])),
        "ql": (([0, .20, .60, 1.0, 1.25, np.nan], [0, .12, .35, .45, .45, np.nan]),
               ([0, .25, .80, 1.4, 2.2, np.nan], [0, .15, .35, .55, .60, np.nan])),
        "w": (([.55, .60, .80, 1.0, 1.1, np.nan], [.10, .15, .30, .30, .35, np.nan]),
              ([.55, .65, 1.4, 2.5, 3.4, np.nan], [.10, .15, .35, .65, .70, np.nan])),
    }
    for ax, (name, label, factor) in zip(axs, conditional):
        group = "thermo" if name in {"thl", "qt", "thv", "ql"} else "default"
        convert = wh_to_z if name == "w" else lambda values: values
        digitized_band(ax, reference, f"f7_{name}")
        # The paper's gray envelopes apply to cloud/core conditional samples.
        # Retain the mean only for the first three scalar panels.
        if name not in {"ql", "w"}:
            ax.plot(factor * convert(d(group, name)), zk, color="tab:blue", label="mean")
        ax.plot(factor * convert(c(group, name)), zk, color="tab:orange", label="cloud")
        ax.plot(factor * convert(co(group, name)), zk, color="tab:green", label="core")
        ax.set_xlabel(label); prof_axis(ax)
    axs[0].legend(fontsize=8)
    fig.suptitle("BOMEX Figure 7 - conditional cloud/core profiles, 3-6 h")
    save(fig, output, "bomex_figure_07_conditional_profiles.png")

    # Figure 8 -- core mass flux and mass-flux flux reconstruction ratio.
    rho = d("thermo", "rho")
    wcore = wh_to_z(co("default", "w"))
    mass_flux = rho * core_area * wcore
    ratios = []
    for scalar, flux in (("qt", d("thermo", "qt_flux")), ("thl", d("thermo", "thl_flux"))):
        core_scalar, mean_scalar = co("thermo", scalar), d("thermo", scalar)
        ratios.append(safe_div(mass_flux * (core_scalar - mean_scalar), wh_to_z(flux)))
    fig, axs = plt.subplots(1, 3, figsize=(10, 4.5), sharey=True)
    digitized_band(axs[0], reference, "f8_mass", 1.1)
    digitized_band(axs[1], reference, "f8_ratio")
    digitized_band(axs[2], reference, "f8_ratio")
    axs[0].plot(mass_flux, zk, color="black"); axs[0].set(xlabel=r"$M_c=\rho a_c w_c$ (kg m$^{-2}$ s$^{-1}$)")
    axs[1].plot(ratios[0], zk, color="black"); axs[1].set(xlabel=r"$M_c(q_{t,c}-\overline{q_t})/w'q_t'$")
    axs[2].plot(ratios[1], zk, color="black"); axs[2].set(xlabel=r"$M_c(\theta_{l,c}-\overline{\theta_l})/w'\theta_l'$")
    for ax in axs: ax.axvline(0, color="0.4", lw=0.8); prof_axis(ax)
    fig.suptitle("BOMEX Figure 8 - core mass-flux diagnostics, 3-6 h")
    save(fig, output, "bomex_figure_08_mass_flux.png")

    # Figure 9 -- diagnosed entrainment/detrainment, neglecting non-conservative sources.
    dmdz = np.gradient(mass_flux, z)
    estimates = []
    for scalar in ("qt", "thl"):
        scalar_core, scalar_env = co("thermo", scalar), d("thermo", scalar)
        transport_gradient = np.gradient(mass_flux * scalar_core, z)
        entrainment = safe_div(transport_gradient - scalar_core * dmdz, scalar_env - scalar_core)
        detrainment = entrainment - dmdz
        estimates.append((safe_div(entrainment, mass_flux), safe_div(detrainment, mass_flux)))
    fig, axs = plt.subplots(1, 2, figsize=(7.4, 4.5), sharey=True)
    for ax, (entrain, detrain), scalar in zip(axs, estimates, ("qt", r"$\theta_l$")):
        ax.plot(1e3 * entrain, zk, label="entrainment")
        ax.plot(1e3 * detrain, zk, label="detrainment")
        ax.axvline(0, color="0.4", lw=0.8); ax.set(xlabel=r"fractional rate (km$^{-1}$)", title=f"from {scalar}")
        prof_axis(ax); ax.legend(fontsize=8)
    fig.suptitle("BOMEX Figure 9 - LES-diagnosed entrainment/detrainment proxy")
    save(fig, output, "bomex_figure_09_entrainment_detrainment.png")

    # Figure 10 -- a simple undiluted buoyancy-integral plume velocity.
    thv_core = co("thermo", "thv")
    buoyancy = G * (thv_core - thv) / thv
    cloud_base = np.argmax(core_area > 1e-4)
    plume_energy = np.zeros_like(z)
    plume_energy[cloud_base:] = 2 * np.maximum.accumulate(np.cumsum(np.maximum(buoyancy[cloud_base:], 0) * np.gradient(z[cloud_base:])))
    w_plume = np.sqrt(plume_energy)
    fig, ax = plt.subplots(figsize=(5, 4.5))
    ax.plot(co("default", "w"), zhk, label="LES cloud core")
    ax.plot(w_plume, zk, "--", label="buoyancy-integral plume")
    ax.set(xlabel="vertical velocity (m s$^{-1}$)"); prof_axis(ax); ax.legend()
    fig.suptitle("BOMEX Figure 10 - core and plume vertical velocity")
    save(fig, output, "bomex_figure_10_plume_velocity.png")

    # Figure 11 -- K = -w'phi'/(d phi/dz).  Use the total LES flux.
    kthl = -safe_div(wh_to_z(d("thermo", "thl_flux")), np.gradient(thl, z))
    kqt = -safe_div(wh_to_z(d("thermo", "qt_flux")), np.gradient(qt, z))
    fig, ax = plt.subplots(figsize=(5, 4.5))
    ax.plot(kthl, zk, label=r"$K_{\theta_l}$")
    ax.plot(kqt, zk, label=r"$K_{q_t}$")
    ax.set(xlabel=r"eddy diffusivity (m$^2$ s$^{-1}$)"); prof_axis(ax); ax.legend()
    fig.suptitle("BOMEX Figure 11 - diagnosed eddy diffusivity, 3-6 h")
    save(fig, output, "bomex_figure_11_eddy_diffusivity.png")

    # Figure 12 -- relative humidity and cloud fraction.  ql_frac is the
    # grid-cell cloud fraction; mask area gives the same cloud diagnostic.
    fig, axs = plt.subplots(1, 2, figsize=(7.4, 4.5), sharey=True)
    axs[0].plot(100 * d("thermo", "rh"), zk); axs[0].set(xlabel="relative humidity (%)")
    axs[1].plot(100 * d("thermo", "ql_frac"), zk, label="ql fraction")
    axs[1].plot(100 * cloud_area, zk, "--", label="ql-mask area")
    axs[1].set(xlabel="cloud fraction (%)"); axs[1].legend(fontsize=8)
    for ax in axs: prof_axis(ax)
    fig.suptitle("BOMEX Figure 12 - relative humidity and cloud fraction, 3-6 h")
    save(fig, output, "bomex_figure_12_rh_cloud_fraction.png")

    # Figure 13 -- instantaneous plan-view cloud occurrence and LWP at the
    # final output time.  A nonzero liquid-water path is the satellite-style
    # top-down definition of cloudy column used in the intercomparison.
    final_time = int(round(default["time"][-1]))
    path_files = sorted(case.glob(f"ql_path.xy.*.{final_time:07d}"))
    figure_count = 12
    if not path_files:
        warnings.warn(
            f"No ql_path xy cross field for t={final_time} s; skipping Figure 13."
        )
    else:
        lwp_raw = np.fromfile(path_files[0], dtype="<f8")
        nxy = int(np.sqrt(lwp_raw.size))
        if nxy * nxy != lwp_raw.size:
            raise ValueError(f"Unexpected LWP field size in {path_files[0].name}")
        lwp = lwp_raw.reshape(nxy, nxy)
        cloud_mask = lwp > 1.e-6
        extent = (0, 12.8, 0, 12.8)
        fig, axs = plt.subplots(1, 2, figsize=(11, 5), sharex=True, sharey=True)
        axs[0].imshow(cloud_mask, origin="lower", extent=extent, cmap="Greys", interpolation="nearest")
        axs[0].set(title="cloud occurrence from above", xlabel="x (km)", ylabel="y (km)")
        visible_lwp = np.ma.masked_less_equal(1e3 * lwp, 1e-3)
        # Clear sky is the dark endpoint of Blues; cloud water brightens through
        # the reversed map, with the largest LWP rendered white.
        axs[1].set_facecolor(plt.get_cmap("Blues")(1.0))
        lwp_image = axs[1].imshow(visible_lwp, origin="lower", extent=extent, cmap="Blues_r",
                                  vmin=0, vmax=np.percentile(1e3 * lwp[cloud_mask], 99.5))
        axs[1].set(title="liquid-water path", xlabel="x (km)")
        fig.colorbar(lwp_image, ax=axs[1], label="LWP (g m$^{-2}$)")
        fig.suptitle(f"BOMEX Figure 13 - instantaneous cloud and LWP, t = {final_time/3600:.1f} h")
        save(fig, output, "bomex_figure_13_final_cloud_lwp.png")
        figure_count += 1

    print(f"Wrote {figure_count} PNGs to {output}")


if __name__ == "__main__":
    main()
