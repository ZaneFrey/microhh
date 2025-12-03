import matplotlib.pyplot as plt
import numpy as np

# ========== DOMAIN INPUTS ==========

# Domain inputs
xsize   = 20 * 1000
ysize   = 20 * 1000
zsize   = 1000
itot    = 256
jtot    = 256
ztot    = 128

# ========== WIND FARM INPUTS ==========

# Number of arrays (wind farms) in the domain
narrays = 2

# Turbine parameters
diameter = 240
height = 150
radius = diameter / 2

# Layout type flags per array (lists of length `narrays`)
flagCenterDom   = [False, False]
flagCenterCoord = [False, False]
flagBottomLeft  = [True, True]
flagStaggeredY  = [False, True]
flagOffset      = [False, False]

# Farm parameters per array
nrows = [4, 5] 
ncols = [4, 5] 
sx    = [8, 7.5] 
sy    = [5, 5]

# If flagBottomLeft: bottom left coordinate of the array 
xstart = [2500, 10000]
ystart = [2500, 10000] 

# If flagCenterCoord: center coordinates of each array
xarrcenter = [2000, 5000]
yarrcenter = [2000, 5000]

# If flagStaggeredY: y-staggering (upward) for every other row of turbines
ystag = [3*diameter, 2*diameter]

# If flagOffset: offset of each row (upward) from the previous
yoffset = [2*diameter, 2*diameter]

# ========== TURBINE POSITIONING ==========

xc_all = []
yc_all = []

for iarr in range(narrays):
    # Per-array settings
    nrows_i = nrows[iarr]
    ncols_i = ncols[iarr]
    sx_i = sx[iarr]
    sy_i = sy[iarr]
    xstart_i = xstart[iarr]
    ystart_i = ystart[iarr]
    flagCenterDom_i = flagCenterDom[iarr]
    flagCenterCoord_i = flagCenterCoord[iarr]
    flagBottomLeft_i = flagBottomLeft[iarr]
    flagStaggeredY_i = flagStaggeredY[iarr]
    flagOffset_i = flagOffset[iarr]
    ystag_i = ystag[iarr]
    yoffset_i = yoffset[iarr]

    # flagBottomLeft defines the bottom-left corner of this array and
    # is mutually exclusive with centering flags for that array.
    if flagBottomLeft_i and (flagCenterDom_i or flagCenterCoord_i):
        raise RuntimeError(
            "flagBottomLeft cannot be True when flagCenterDom or "
            "flagCenterCoord is True for array {}".format(iarr)
        )

    dx_turb = sx_i * diameter
    dy_turb = sy_i * diameter

    x_positions = xstart_i + np.arange(nrows_i) * dx_turb
    y_positions = ystart_i + np.arange(ncols_i) * dy_turb

    xc_grid, yc_grid = np.meshgrid(x_positions, y_positions)

    # Staggering/offset cannot both be enabled for a single array.
    if flagStaggeredY_i and flagOffset_i:
        raise RuntimeError(
            "flagStaggeredY and flagOffset cannot both be True for array {}".format(
                iarr
            )
        )

    # Apply stagger/offset in y based on streamwise row index (x-direction)
    if flagStaggeredY_i:
        row_indices_x = np.arange(xc_grid.shape[1])
        y_shift = (row_indices_x % 2) * ystag_i
        yc_grid = yc_grid + y_shift
    elif flagOffset_i:
        row_indices_x = np.arange(xc_grid.shape[1])
        y_shift = row_indices_x * yoffset_i
        yc_grid = yc_grid + y_shift

    # After any staggering/offset, center the final turbine array
    # according to the selected centering flag (if enabled).
    if (flagCenterDom_i or flagCenterCoord_i) and not flagBottomLeft_i:
        x_min = xc_grid.min()
        x_max = xc_grid.max()
        y_min = yc_grid.min()
        y_max = yc_grid.max()
        x_center_array = 0.5 * (x_min + x_max)
        y_center_array = 0.5 * (y_min + y_max)

        if flagCenterDom_i:
            x_target = 0.5 * xsize
            y_target = 0.5 * ysize
        else:
            x_target = xarrcenter[iarr]
            y_target = yarrcenter[iarr]

        dx_center = x_target - x_center_array
        dy_center = y_target - y_center_array
        xc_grid = xc_grid + dx_center
        yc_grid = yc_grid + dy_center

    xc_all.append(xc_grid.ravel())
    yc_all.append(yc_grid.ravel())

if narrays == 1:
    xc = xc_all[0]
    yc = yc_all[0]
else:
    xc = np.concatenate(xc_all)
    yc = np.concatenate(yc_all)

# Print grid resolution diagnostic and turbine locations in .ini-style format
print()
dy_grid = float(ysize) / jtot
ny_per_diam = diameter / dy_grid
print("Grid points in y across rotor = {:.2f}".format(ny_per_diam))
print("xc=" + ",".join("{:.0f}".format(x) for x in xc))
print("yc=" + ",".join("{:.0f}".format(y) for y in yc))

# ========== PLOTTING ==========

fig, ax = plt.subplots()
ax.set_aspect("equal", adjustable="box")

# Plot domain extent as a black rectangle based on xsize and ysize
ax.plot(
    [0, xsize, xsize, 0, 0],
    [0, 0, ysize, ysize, 0],
    color="k",
    linewidth=1.5,
)

# Plot the numerical grid based on itot and jtot as faint gray lines
dx = float(xsize) / itot
dy = float(ysize) / jtot

for i in range(itot + 1):
    x = i * dx
    ax.plot([x, x], [0, ysize], color="0.8", linewidth=0.5, zorder=0)

for j in range(jtot + 1):
    y = j * dy
    ax.plot([0, xsize], [y, y], color="0.8", linewidth=0.5, zorder=0)

# Extract the turbine center locations and plot center locations as small red dot
for x_t, y_t in zip(xc, yc):
    # Turbine center
    ax.plot(x_t, y_t, "ro", markersize=1, zorder=2)

    # Plot extent of disk
    ax.plot(
        [x_t, x_t],
        [y_t - radius, y_t + radius],
        color="r",
        linewidth=1.0,
        zorder=2,
    )

    # Plot circle around turbine swept extent
    circle = plt.Circle(
        (x_t, y_t),
        radius,
        edgecolor="r",
        facecolor="none",
        linestyle="--",
        linewidth=0.5,
        alpha=0.3,
        zorder=1,
    )
    ax.add_patch(circle)

ax.set_xlim(0, xsize)
ax.set_ylim(0, ysize)
ax.set_xlabel("x [m]")
ax.set_ylabel("y [m]")

plt.tight_layout()
plt.show()
