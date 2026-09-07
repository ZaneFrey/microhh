"""Create the neutral Wu and Porté-Agel (2011) initial velocity profiles."""

from configparser import ConfigParser
from pathlib import Path

import netCDF4 as nc
import numpy as np


case_name = Path(__file__).stem.removesuffix("_input")
case_dir = Path(__file__).resolve().parent

ini = ConfigParser()
ini.read(case_dir / f"{case_name}.ini")

ktot = ini.getint("grid", "ktot")
zsize = ini.getfloat("grid", "zsize")
z0m = ini.getfloat("boundary", "z0m")
hub_height = ini.getfloat("windfarm", "hubheight")

# Wu and Porté-Agel (2011), Sect. 3: neutral boundary layer with u_* =
# 0.102 m s-1, z0 = 0.03 mm, U_infinity = 2.8 m s-1, and U_hub = 2.2 m s-1.
# The logarithmic profile is used only to initialise the pressure-driven LES;
# the specified pressure gradient and surface roughness determine its evolved
# statistically stationary state.  It is scaled to the reported hub-height
# velocity and smoothly blended to the reported free-stream velocity at the
# boundary-layer top.
ustar = 0.102
kappa = 0.4
freestream_velocity = 2.8
hub_velocity = 2.2

dz = zsize / ktot
z = np.linspace(0.5 * dz, zsize - 0.5 * dz, ktot)
u_log = ustar / kappa * np.log(z / z0m)
u_hub_log = ustar / kappa * np.log(hub_height / z0m)
u = u_log * (hub_velocity / u_hub_log)

# Cubic Hermite blending keeps the profile and its slope continuous at the hub
# and at the boundary-layer top.
blend = np.clip((z - hub_height) / (zsize - hub_height), 0.0, 1.0)
blend = blend * blend * (3.0 - 2.0 * blend)
u = (1.0 - blend) * u + blend * freestream_velocity
v = np.zeros_like(u)

with nc.Dataset(case_dir / f"{case_name}_input.nc", "w", format="NETCDF4") as dataset:
    dataset.createDimension("z", ktot)
    dataset.createVariable("z", "f8", ("z",))[:] = z

    initial = dataset.createGroup("init")
    initial.createVariable("u", "f8", ("z",))[:] = u
    initial.createVariable("v", "f8", ("z",))[:] = v
