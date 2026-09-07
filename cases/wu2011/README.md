# Wu and Porté-Agel (2011) wind-turbine wake case

## Case description

This case is a dimensional MicroHH realization of the neutral wind-turbine wake experiment described by Wu and Porté-Agel (2011), *Large-Eddy Simulation of Wind-Turbine Wakes: Evaluation of Turbine Parametrisations*. It is intended to validate the actuator-disk model (ADM) and actuator-disk model with rotation (ADM-R).

The pressure-driven boundary layer has a depth of `0.46 m`, friction velocity `u_* = 0.102 m s-1`, and aerodynamic roughness length `z0 = 0.03 mm`. The initial velocity profile is scaled to the reported hub-height velocity of `2.2 m s-1` and blended toward the reported free-stream velocity of `2.8 m s-1` near the top of the domain.

The domain is `8.64 x 0.72 x 0.46 m^3` on a `576 x 48 x 64` grid, giving horizontal grid spacings of `0.015 m`. A single turbine with diameter `D = 0.15 m` and hub height `zh = 0.125 m` is located at `(x, y) = (0.9, 0.36) m`. The default `swwf=admr` configuration uses `C_T' = 4/3`, `C_P' = 0.8`, and a tip-speed ratio of `4`; set `swwf=adm` to run the non-rotating actuator-disk model.

The simulation runs for `600 s`. Flow statistics, three-dimensional velocity dumps, and hub-height and rotor-plane cross-sections are sampled every `10 s`; turbine data are sampled every time step, and restart files are saved every `60 s`.