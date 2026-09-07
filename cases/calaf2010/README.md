# Fully developed wind-turbine array boundary layer

## Case description

This case is a dimensional MicroHH realization of the baseline A1 wind-turbine-array boundary layer from *Calaf, Meneveau, and Meyers (2010): Large eddy simulation study of fully developed wind-turbine array boundary layers, Physics of Fluids 22, 015110*. It uses a neutral, doubly periodic domain with `H = 1000 m`, `Lx = Ly = pi H`, a `128^3` grid, and 24 actuator disks arranged in four streamwise columns and six spanwise rows. The turbine diameter and hub height are both `100 m`, the lower-surface roughness is `z0 = 0.1 m`, and `C_T' = 4/3`.

## Running the case

Generate the input NetCDF file and launch the case with a MicroHH executable containing the WindFarm dynamic-yaw implementation:

```bash
python3 calaf2010_input.py
../../build-develop_windfarm/microhh init calaf2010
../../build-develop_windfarm/microhh run calaf2010
```
