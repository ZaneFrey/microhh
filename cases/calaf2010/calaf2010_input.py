"""Create initial and geostrophic profiles for the Calaf et al. (2010) case."""

from configparser import ConfigParser
from pathlib import Path

import netCDF4 as nc
import numpy as np


case_name = Path(__file__).stem.removesuffix("_input")
case_dir = Path(__file__).resolve().parent
ini_path = case_dir / f"{case_name}.ini"

ini = ConfigParser()
if not ini.read(ini_path):
    raise FileNotFoundError(f"Could not read {ini_path}")

ktot = ini.getint("grid", "ktot")
zsize = ini.getfloat("grid", "zsize")
z0m = ini.getfloat("boundary", "z0m")
forcing = ini.get("force", "swlspres")
fc = ini.getfloat("force", "fc")

if ktot <= 0 or zsize <= 0.0 or z0m <= 0.0:
    raise ValueError("ktot, zsize, and z0m must be positive")
if forcing != "geo" or fc <= 0.0:
    raise ValueError("This input generator expects positive-f Coriolis/geostrophic forcing")

# The Appendix uses f=9.34e-5 s-1 at 40 degrees latitude and discusses a
# geostrophic vector (U_G, V_G)=(10, -5) m s-1 for the no-farm reference with
# u_*=0.45 m s-1.  The negative transverse component makes x the direction of
# the near-surface wind, as assumed in the Appendix.  Once turbines are active,
# u_* and the boundary-layer depth become outcomes rather than fixed inputs.
ustar_reference = 0.45
u_geo_value = 10.0
v_geo_value = -5.0
kappa = 0.4

dz = zsize / ktot
z = (np.arange(ktot, dtype=np.float64) + 0.5) * dz

# Initialize with the neutral, no-turbine logarithmic equilibrium profile from
# Eq. (5). The actuator disks and geostrophic forcing then spin this profile up to
# the fully developed wind-turbine-array state. Random horizontal perturbations
# configured in the INI seed the resolved turbulence during `microhh init`.
u = ustar_reference / kappa * np.log(z / z0m)
v = np.zeros_like(u)
u_geo = np.full_like(z, u_geo_value)
v_geo = np.full_like(z, v_geo_value)

output_path = case_dir / f"{case_name}_input.nc"
with nc.Dataset(output_path, "w", format="NETCDF4") as dataset:
    dataset.createDimension("z", ktot)

    z_var = dataset.createVariable("z", "f8", ("z",))
    z_var.units = "m"
    z_var.long_name = "height at cell center"
    z_var[:] = z

    initial = dataset.createGroup("init")

    u_var = initial.createVariable("u", "f8", ("z",))
    u_var.units = "m s-1"
    u_var.long_name = "initial streamwise velocity"
    u_var[:] = u

    v_var = initial.createVariable("v", "f8", ("z",))
    v_var.units = "m s-1"
    v_var.long_name = "initial spanwise velocity"
    v_var[:] = v

    u_geo_var = initial.createVariable("u_geo", "f8", ("z",))
    u_geo_var.units = "m s-1"
    u_geo_var.long_name = "streamwise geostrophic wind"
    u_geo_var[:] = u_geo

    v_geo_var = initial.createVariable("v_geo", "f8", ("z",))
    v_geo_var.units = "m s-1"
    v_geo_var.long_name = "spanwise geostrophic wind"
    v_geo_var[:] = v_geo

print(
    f"Wrote {output_path.name}: {ktot} levels, dz={dz:.6g} m, "
    f"reference u_*={ustar_reference:.6g} m s-1, "
    f"(U_G, V_G)=({u_geo_value:.6g}, {v_geo_value:.6g}) m s-1"
)
