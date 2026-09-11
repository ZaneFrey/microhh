#!/usr/bin/env python3
"""Compare five BOMEX-WF advection configurations.

The script reproduces comparison equivalents of Siebesma et al. (2003)
Figures 1-8 and 13. Figures 9-12 are intentionally omitted. Each diagnostic
panel contains one curve per advection configuration. Figure 13 contains one
cloud-occurrence/LWP row per configuration.
"""

from __future__ import annotations

import argparse
import configparser
from collections import OrderedDict
from pathlib import Path
import warnings

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from analyze_bomex_intercomparison import (
    AVG_END,
    AVG_START,
    CP,
    KM,
    LV,
    RHO_REF,
    field,
    load,
    load_restart_series,
    mean_last3,
    reliable_flux,
    safe_div,
)


CASE_SPECS = OrderedDict([
    ("swadvec2", ("2", Path("old/swadvec2"))),
    ("swadvec2i5", ("2i5", Path("old/swadvec2i5"))),
    ("swadvec2i5limited", ("2i5 + limiter", Path("old/swadvec2i5lim"))),
    ("swadvec2i62", ("2i62", Path("old/swadvec2i62"))),
    ("swadvec2i62limited", ("2i62 + limiter", Path("old/swadvec2i62lim"))),
])


def style() -> None:
    plt.rcParams.update({
        "figure.dpi": 140,
        "savefig.dpi": 180,
        "font.size": 10,
        "axes.grid": True,
        "grid.alpha": 0.25,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "lines.linewidth": 1.8,
    })


def save(fig: plt.Figure, output: Path, name: str) -> None:
    fig.tight_layout()
    fig.savefig(output / name, bbox_inches="tight")
    plt.close(fig)


def prof_axis(ax: plt.Axes, ylabel: bool = True) -> None:
    ax.set_ylim(0, 3.5)
    ax.set_ylabel("height (km)" if ylabel else "")


def case_legend(fig: plt.Figure, axes, extra_handles=()) -> None:
    handles, labels = axes.flat[0].get_legend_handles_labels()
    for handle, label in extra_handles:
        handles.append(handle)
        labels.append(label)
    fig.legend(handles, labels, loc="upper center", ncol=3,
               bbox_to_anchor=(0.5, 0.995), frameon=False)


def center_w(values: np.ndarray) -> np.ndarray:
    return 0.5 * (values[..., :-1] + values[..., 1:])


def load_run(path: Path) -> dict:
    run = {
        "path": path,
        "default": load_restart_series(path, "bomex-wf.default"),
        "cloud": load_restart_series(path, "bomex-wf.ql"),
        "core": load_restart_series(path, "bomex-wf.qlcore"),
        "initial": load(path / "bomex-wf_input.nc"),
    }
    return run


def averaged(run: dict, sample: str, group: str, name: str) -> np.ndarray:
    return mean_last3(run[sample], group, name)


def load_lwp(run: dict, nx: int, ny: int) -> tuple[np.ndarray, int]:
    files = sorted(
        run["path"].glob("ql_path.xy.*.[0-9][0-9][0-9][0-9][0-9][0-9][0-9]"),
        key=lambda path: int(path.name.rsplit(".", 1)[1]),
    )
    if not files:
        raise FileNotFoundError(f"No ql_path xy fields in {run['path']}")

    path = files[-1]
    count = nx * ny
    size = path.stat().st_size
    if size == 8 * count:
        dtype = "<f8"
    elif size == 4 * count:
        dtype = "<f4"
    else:
        raise ValueError(f"Unexpected LWP field size for {path}: {size} bytes")

    return np.fromfile(path, dtype=dtype).reshape(ny, nx), int(path.name.rsplit(".", 1)[1])


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case-dir", type=Path, default=Path(__file__).resolve().parent)
    parser.add_argument("--output-dir", type=Path)
    args = parser.parse_args()

    case_dir = args.case_dir.resolve()
    output = (args.output_dir.resolve() if args.output_dir
              else case_dir / "outputs" / "advectest")
    output.mkdir(parents=True, exist_ok=True)

    style()
    colors = plt.get_cmap("tab10").colors[:len(CASE_SPECS)]
    runs = OrderedDict()
    for color, (name, (label, relative_path)) in zip(colors, CASE_SPECS.items()):
        path = case_dir / relative_path
        if not path.is_dir():
            raise FileNotFoundError(f"Missing advection run directory: {path}")
        run = load_run(path)
        run.update(label=label, color=color)
        runs[name] = run

    z_reference = next(iter(runs.values()))["default"]["z"]
    zh_reference = next(iter(runs.values()))["default"]["zh"]
    for name, run in runs.items():
        if not np.allclose(run["default"]["z"], z_reference):
            raise ValueError(f"Cell-center grid differs for {name}")
        if not np.allclose(run["default"]["zh"], zh_reference):
            raise ValueError(f"Cell-edge grid differs for {name}")
    zk, zhk = z_reference * KM, zh_reference * KM

    # Figure 1: identical initial profiles are retained as a reproducibility check.
    fig, axs = plt.subplots(1, 4, figsize=(12, 4.5), sharey=True)
    initial_panels = [
        ("thl", 1.0, r"$\theta_l$ (K)"),
        ("qt", 1e3, r"$q_t$ (g kg$^{-1}$)"),
        ("u", 1.0, "u (m s$^{-1}$)"),
        ("v", 1.0, "v (m s$^{-1}$)"),
    ]
    for run in runs.values():
        for ax, (name, factor, xlabel) in zip(axs, initial_panels):
            ax.plot(factor * field(run["initial"], "init", name), zk,
                    color=run["color"], label=run["label"])
            ax.set_xlabel(xlabel)
    for index, ax in enumerate(axs):
        prof_axis(ax, index == 0)
    fig.suptitle("BOMEX advection comparison - Figure 1 initial profiles", y=1.06)
    case_legend(fig, axs)
    save(fig, output, "advec_figure_01_initial_profiles.png")

    # Figure 2: ql_cover is total projected cloud cover, not ql_frac(z).
    fig, axs = plt.subplots(1, 3, figsize=(12, 4.2))
    for run in runs.values():
        default = run["default"]
        time_h = default["time"] / 3600.0
        tke = 0.5 * (
            field(default, "default", "u_2")
            + field(default, "default", "v_2")
            + center_w(field(default, "default", "w_2"))
        )
        axs[0].plot(time_h, 100 * field(default, "thermo", "ql_cover"),
                    color=run["color"], label=run["label"])
        axs[1].plot(time_h, 1e3 * field(default, "thermo", "ql_path"),
                    color=run["color"], label=run["label"])
        axs[2].plot(time_h, np.trapezoid(tke, z_reference, axis=1),
                    color=run["color"], label=run["label"])
    ylabels = ["cloud cover (%)", "LWP (g m$^{-2}$)",
               r"$\int$ TKE dz (m$^3$ s$^{-2}$)"]
    for ax, ylabel in zip(axs, ylabels):
        ax.axvspan(3, 6, color="0.9", zorder=0)
        ax.set(xlim=(0, 6), xlabel="time (h)", ylabel=ylabel)
    fig.suptitle("BOMEX advection comparison - Figure 2 domain evolution", y=1.06)
    case_legend(fig, axs)
    save(fig, output, "advec_figure_02_time_series.png")

    # Figure 3: mean profiles over hours 3-6 plus the common initial profile.
    fig, axs = plt.subplots(1, 5, figsize=(15, 4.5), sharey=True)
    panels = [
        ("thermo", "thl", 1.0, r"$\theta_l$ (K)", "thl"),
        ("thermo", "qt", 1e3, r"$q_t$ (g kg$^{-1}$)", "qt"),
        ("default", "u", 1.0, "u (m s$^{-1}$)", "u"),
        ("default", "v", 1.0, "v (m s$^{-1}$)", "v"),
        ("thermo", "ql", 1e3, r"$q_l$ (g kg$^{-1}$)", None),
    ]
    for run in runs.values():
        for ax, (group, name, factor, xlabel, _) in zip(axs, panels):
            ax.plot(factor * averaged(run, "default", group, name), zk,
                    color=run["color"], label=run["label"])
            ax.set_xlabel(xlabel)
    baseline_initial = next(iter(runs.values()))["initial"]
    initial_handle = None
    for ax, (_, _, factor, _, initial_name) in zip(axs, panels):
        initial_values = (factor * field(baseline_initial, "init", initial_name)
                          if initial_name else np.zeros_like(z_reference))
        handle, = ax.plot(initial_values, zk, "k--", lw=1.2)
        initial_handle = handle
    for index, ax in enumerate(axs):
        prof_axis(ax, index == 0)
    fig.suptitle("BOMEX advection comparison - Figure 3 mean profiles, 3-6 h", y=1.06)
    case_legend(fig, axs, ((initial_handle, "initial"),))
    save(fig, output, "advec_figure_03_mean_profiles.png")

    # Figure 4: scalar fluxes are converted to the paper's energetic units.
    fig, axs = plt.subplots(1, 5, figsize=(15, 4.5), sharey=True)
    for run in runs.values():
        ql_flux, _ = reliable_flux(
            averaged(run, "default", "thermo", "ql_flux"),
            averaged(run, "default", "thermo", "ql_w"), "ql", 1e-2)
        thv_flux, _ = reliable_flux(
            averaged(run, "default", "thermo", "thv_flux"),
            averaged(run, "default", "thermo", "thv_w"), "thv", 10.0)
        profiles = [
            RHO_REF * LV * averaged(run, "default", "thermo", "qt_flux"),
            RHO_REF * CP * averaged(run, "default", "thermo", "thl_flux"),
            RHO_REF * LV * ql_flux,
            RHO_REF * CP * thv_flux,
            averaged(run, "default", "default", "u_flux"),
        ]
        for ax, values in zip(axs, profiles):
            ax.plot(values, zhk, color=run["color"], label=run["label"])
    xlabels = [
        r"$\rho L_v\overline{w'q_t'}$ (W m$^{-2}$)",
        r"$\rho c_p\overline{w'\theta_l'}$ (W m$^{-2}$)",
        r"$\rho L_v\overline{w'q_l'}$ (W m$^{-2}$)",
        r"$\rho c_p\overline{w'\theta_v'}$ (W m$^{-2}$)",
        r"$\overline{u'w'}$ (m$^2$ s$^{-2}$)",
    ]
    for index, (ax, xlabel) in enumerate(zip(axs, xlabels)):
        ax.set_xlabel(xlabel)
        prof_axis(ax, index == 0)
    fig.suptitle("BOMEX advection comparison - Figure 4 turbulent fluxes, 3-6 h", y=1.06)
    case_legend(fig, axs)
    save(fig, output, "advec_figure_04_flux_profiles.png")

    # Figure 5: TKE and vertical-velocity variance.
    fig, axs = plt.subplots(1, 2, figsize=(8.5, 4.5), sharey=True)
    for run in runs.values():
        default = run["default"]
        use = (default["time"] >= AVG_START) & (default["time"] <= AVG_END)
        tke = 0.5 * (
            field(default, "default", "u_2")
            + field(default, "default", "v_2")
            + center_w(field(default, "default", "w_2"))
        )
        axs[0].plot(np.nanmean(tke[use], axis=0), zk,
                    color=run["color"], label=run["label"])
        axs[1].plot(averaged(run, "default", "default", "w_2"), zhk,
                    color=run["color"], label=run["label"])
    axs[0].set_xlabel("TKE (m$^2$ s$^{-2}$)")
    axs[1].set_xlabel(r"$\overline{w'^2}$ (m$^2$ s$^{-2}$)")
    prof_axis(axs[0], True)
    prof_axis(axs[1], False)
    fig.suptitle("BOMEX advection comparison - Figure 5 turbulence, 3-6 h", y=1.06)
    case_legend(fig, axs)
    save(fig, output, "advec_figure_05_tke_wvariance.png")

    # Figure 6: separate panels keep five configuration curves per diagnostic.
    fig, axs = plt.subplots(1, 2, figsize=(8.5, 4.5), sharey=True)
    for run in runs.values():
        cloud_area = averaged(run, "cloud", "default", "area")
        core_area = averaged(run, "core", "default", "area")
        axs[0].plot(100 * cloud_area, zk, color=run["color"], label=run["label"])
        axs[1].plot(100 * core_area, zk, color=run["color"], label=run["label"])
    axs[0].set(xlabel="cloud fractional area (%)", title=r"cloud: $q_l>0$")
    axs[1].set(xlabel="core fractional area (%)", title=r"core: $q_l>0$, $b'>0$")
    prof_axis(axs[0], True)
    prof_axis(axs[1], False)
    fig.suptitle("BOMEX advection comparison - Figure 6 cloud and core cover, 3-6 h", y=1.06)
    case_legend(fig, axs)
    save(fig, output, "advec_figure_06_cloud_core_cover.png")

    # Figure 7: cloud and core conditional profiles occupy separate rows.
    fig, axs = plt.subplots(2, 5, figsize=(15, 8), sharey=True)
    conditional = [
        ("thermo", "thl", 1.0, r"$\theta_l$ (K)"),
        ("thermo", "qt", 1e3, r"$q_t$ (g kg$^{-1}$)"),
        ("thermo", "thv", 1.0, r"$\theta_v$ (K)"),
        ("thermo", "ql", 1e3, r"$q_l$ (g kg$^{-1}$)"),
        ("default", "w", 1.0, "w (m s$^{-1}$)"),
    ]
    for run in runs.values():
        for row, sample in enumerate(("cloud", "core")):
            for column, (group, name, factor, xlabel) in enumerate(conditional):
                values = averaged(run, sample, group, name)
                if name == "w":
                    values = center_w(values)
                axs[row, column].plot(factor * values, zk,
                                      color=run["color"], label=run["label"])
                axs[row, column].set_xlabel(xlabel)
    axs[0, 0].set_title("cloud conditional")
    axs[1, 0].set_title("core conditional")
    for row in range(2):
        for column in range(5):
            prof_axis(axs[row, column], column == 0)
    fig.suptitle("BOMEX advection comparison - Figure 7 conditional profiles, 3-6 h", y=1.03)
    case_legend(fig, axs)
    save(fig, output, "advec_figure_07_conditional_profiles.png")

    # Figure 8: core mass flux and reconstruction ratios.
    fig, axs = plt.subplots(1, 3, figsize=(11, 4.5), sharey=True)
    for run in runs.values():
        core_area = averaged(run, "core", "default", "area")
        rho = averaged(run, "default", "thermo", "rho")
        wcore = center_w(averaged(run, "core", "default", "w"))
        mass_flux = rho * core_area * wcore
        ratios = []
        for scalar in ("qt", "thl"):
            scalar_core = averaged(run, "core", "thermo", scalar)
            scalar_mean = averaged(run, "default", "thermo", scalar)
            total_flux = center_w(averaged(run, "default", "thermo", f"{scalar}_flux"))
            ratios.append(safe_div(mass_flux * (scalar_core - scalar_mean), total_flux))
        for ax, values in zip(axs, (mass_flux, *ratios)):
            ax.plot(values, zk, color=run["color"], label=run["label"])
    xlabels = [
        r"$M_c=\rho a_cw_c$ (kg m$^{-2}$ s$^{-1}$)",
        r"$M_c(q_{t,c}-\overline{q_t})/\overline{w'q_t'}$",
        r"$M_c(\theta_{l,c}-\overline{\theta_l})/\overline{w'\theta_l'}$",
    ]
    for index, (ax, xlabel) in enumerate(zip(axs, xlabels)):
        ax.axvline(0, color="0.5", lw=0.8)
        ax.set_xlabel(xlabel)
        prof_axis(ax, index == 0)
    fig.suptitle("BOMEX advection comparison - Figure 8 core mass-flux diagnostics, 3-6 h", y=1.06)
    case_legend(fig, axs)
    save(fig, output, "advec_figure_08_mass_flux.png")

    # Figure 13: final projected cloud mask and LWP for all configurations.
    config = configparser.ConfigParser()
    config.read(case_dir / "bomex-wf.ini")
    nx = config.getint("grid", "itot")
    ny = config.getint("grid", "jtot")
    extent = (0, config.getfloat("grid", "xsize") * KM,
              0, config.getfloat("grid", "ysize") * KM)
    lwp_fields = OrderedDict()
    positive_values = []
    for name, run in runs.items():
        lwp, time = load_lwp(run, nx, ny)
        lwp_fields[name] = (1e3 * lwp, time)
        positive_values.append((1e3 * lwp)[lwp > 0])
    positive_values = np.concatenate([values for values in positive_values if values.size])
    vmax = np.percentile(positive_values, 99.5)

    fig, axs = plt.subplots(len(runs), 2, figsize=(11, 22), sharex=True, sharey=True,
                            layout="constrained")
    lwp_image = None
    for row, (name, run) in enumerate(runs.items()):
        lwp, time = lwp_fields[name]
        cloud_mask = lwp > 1e-3
        axs[row, 0].imshow(cloud_mask, origin="lower", extent=extent, cmap="Greys",
                          interpolation="nearest", vmin=0, vmax=1)
        visible_lwp = np.ma.masked_less_equal(lwp, 1e-3)
        axs[row, 1].set_facecolor(plt.get_cmap("Blues")(1.0))
        lwp_image = axs[row, 1].imshow(
            visible_lwp, origin="lower", extent=extent, cmap="Blues_r",
            interpolation="nearest", vmin=0, vmax=vmax)
        axs[row, 0].set_ylabel(f"{run['label']}\ny (km)")
        axs[row, 0].set_title(f"cloud occurrence, t={time/3600:.1f} h")
        axs[row, 1].set_title(f"LWP, t={time/3600:.1f} h")
    axs[-1, 0].set_xlabel("x (km)")
    axs[-1, 1].set_xlabel("x (km)")
    fig.colorbar(lwp_image, ax=axs[:, 1], label="LWP (g m$^{-2}$)", shrink=0.8)
    fig.suptitle("BOMEX advection comparison - Figure 13 final cloud fields")
    fig.savefig(output / "advec_figure_13_final_cloud_lwp.png", bbox_inches="tight")
    plt.close(fig)

    print("\nHours 3-6 summary")
    print(f"{'configuration':24s} {'cover (%)':>12s} {'LWP (g m-2)':>14s}")
    for run in runs.values():
        default = run["default"]
        use = (default["time"] >= AVG_START) & (default["time"] <= AVG_END)
        cover = 100 * np.mean(field(default, "thermo", "ql_cover")[use])
        lwp = 1e3 * np.mean(field(default, "thermo", "ql_path")[use])
        print(f"{run['label']:24s} {cover:12.4f} {lwp:14.4f}")
    print(f"\nWrote 9 comparison PNGs to {output}")


if __name__ == "__main__":
    main()
