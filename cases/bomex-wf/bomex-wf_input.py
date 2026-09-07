import numpy as np
import netCDF4 as nc

float_type = "f8"

# Get number of vertical levels and size from .ini file
with open('bomex-wf.ini') as f:
    for line in f:
        if(line.split('=')[0]=='ktot'):
            kmax = int(line.split('=')[1])
        if(line.split('=')[0]=='zsize'):
            zsize = float(line.split('=')[1])


# Discretize grid with hyperbolic tangent function
def target_dz(height):
    dz_min = 5.0
    dz_max = 25.0
    transition_height = 1588.0
    transition_width = 500.0

    return (
        dz_min
        + 0.5 * (dz_max - dz_min)
        * (1.0 + np.tanh(
            (height - transition_height) / transition_width
        ))
    )

# Integrate dz/dk = target_dz(z) using RK4.
zh = np.zeros(kmax + 1)

for k in range(kmax):
    height = zh[k]

    rk1 = target_dz(height)
    rk2 = target_dz(height + 0.5 * rk1)
    rk3 = target_dz(height + 0.5 * rk2)
    rk4 = target_dz(height + rk3)

    zh[k + 1] = height + (
        rk1 + 2.0 * rk2 + 2.0 * rk3 + rk4
    ) / 6.0

# The unscaled grid ends at 3998.77 m. Apply a negligible
# 0.0307% correction so the upper boundary is exactly zsize.
zh *= zsize / zh[-1]

# MicroHH expects cell-center heights in the input file.
z = 0.5 * (zh[:-1] + zh[1:])
dz = np.diff(zh)


# set the height
thl   = np.zeros(np.size(z))
qt    = np.zeros(np.size(z))
u     = np.zeros(np.size(z))
v     = np.zeros(np.size(z))
ugeo  = np.zeros(np.size(z))
vgeo  = np.zeros(np.size(z))
wls   = np.zeros(np.size(z))
thlls = np.zeros(np.size(z))
qtls  = np.zeros(np.size(z))

for k in range(kmax):
    zk = z[k]

    # Liquid-water potential temperature
    if zk <= 520.0:
        thl[k] = 298.7
    elif zk <= 1480.0:
        thl[k] = (
            298.7
            + (zk - 520.0) * (302.4 - 298.7) / (1480.0 - 520.0)
        )
    elif zk <= 2000.0:
        thl[k] = (
            302.4
            + (zk - 1480.0) * (308.2 - 302.4) / (2000.0 - 1480.0)
        )
    elif zk <= 3000.0:
        thl[k] = (
            308.2
            + (zk - 2000.0) * (311.85 - 308.2) / (3000.0 - 2000.0)
        )
    else:
        # Continue the 3.65 K/km free-tropospheric gradient.
        thl[k] = 311.85 + (zk - 3000.0) * 3.65e-3

    # Total specific humidity
    if zk <= 520.0:
        qt[k] = 1.0e-3 * (
            17.0 + zk * (16.3 - 17.0) / 520.0
        )
    elif zk <= 1480.0:
        qt[k] = 1.0e-3 * (
            16.3
            + (zk - 520.0) * (10.7 - 16.3) / (1480.0 - 520.0)
        )
    elif zk <= 2000.0:
        qt[k] = 1.0e-3 * (
            10.7
            + (zk - 1480.0) * (4.2 - 10.7) / (2000.0 - 1480.0)
        )
    elif zk <= 3000.0:
        qt[k] = 1.0e-3 * (
            4.2
            + (zk - 2000.0) * (3.0 - 4.2) / (3000.0 - 2000.0)
        )
    else:
        # Continue the -1.2 g/kg/km gradient.
        qt[k] = 1.0e-3 * (
            3.0 + (zk - 3000.0) * (-1.2) / 1000.0
        )

    # Initial horizontal wind
    if zk <= 700.0:
        u[k] = -8.75
    else:
        # Continue the original 1.8e-3 s^-1 shear.
        u[k] = -8.75 + (zk - 700.0) * 1.8e-3

    # Geostrophic wind
    ugeo[k] = -10.0 + 1.8e-3 * zk

    # Large-scale subsidence: zero above 2100 m
    if zk <= 1500.0:
        wls[k] = zk * (-0.65) / 1500.0
    elif zk <= 2100.0:
        wls[k] = -0.65 + (zk - 1500.0) * 0.65 / 600.0
    else:
        wls[k] = 0.0

    # Large-scale temperature tendency: zero above 3000 m
    if zk <= 1500.0:
        thlls[k] = -2.0
    elif zk <= 3000.0:
        thlls[k] = -2.0 + (zk - 1500.0) * 2.0 / 1500.0
    else:
        thlls[k] = 0.0

    # Large-scale moisture tendency: zero above 500 m
    if zk <= 300.0:
        qtls[k] = -1.2
    elif zk <= 500.0:
        qtls[k] = -1.2 + (zk - 300.0) * 1.2 / 200.0
    else:
        qtls[k] = 0.0
        

# normalize profiles to SI
#qtls /= 1000.  # from g/kg to kg/kg
wls  /= 100.   # from cm/s to m/s
thlls  /= 86400. # from K/d to K/s
qtls *= 1.e-8

nc_file = nc.Dataset("bomex-wf_input.nc", mode="w", datamodel="NETCDF4", clobber=True)
nc_file.createDimension("z", kmax)
nc_z = nc_file.createVariable("z", float_type, ("z"))

nc_group_init = nc_file.createGroup("init");
nc_thl   = nc_group_init.createVariable("thl"   , float_type, ("z"))
nc_qt    = nc_group_init.createVariable("qt"    , float_type, ("z"))
nc_u     = nc_group_init.createVariable("u"     , float_type, ("z"))
nc_v     = nc_group_init.createVariable("v"     , float_type, ("z"))
nc_ugeo  = nc_group_init.createVariable("u_geo" , float_type, ("z"))
nc_vgeo  = nc_group_init.createVariable("v_geo" , float_type, ("z"))
nc_wls   = nc_group_init.createVariable("w_ls"  , float_type, ("z"))
nc_thlls = nc_group_init.createVariable("thl_ls", float_type, ("z"))
nc_qtls  = nc_group_init.createVariable("qt_ls" , float_type, ("z"))

nc_z    [:] = z    [:]
nc_thl  [:] = thl  [:]
nc_qt   [:] = qt   [:]
nc_u    [:] = u    [:]
nc_ugeo [:] = ugeo [:]
nc_v    [:] = v    [:]
nc_vgeo [:] = vgeo [:]
nc_wls  [:] = wls  [:]
nc_thlls[:] = thlls[:]
nc_qtls [:] = qtls [:]

nc_file.close()
