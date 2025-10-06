# CM1 Heterogeneous Initialization

Modify CM1 v21.1 to initialize from a **3-D heterogeneous analysis** (e.g., HRRR analysis). This adds a new initialization path when you set `isnd=67` that:
- reads a NetCDF file with 3-D fields on a source grid named "input_sounding.nc",
- horizontally resamples each source level to the CM1 horizontal grid specifications,
- vertically and horizontally interpolates to CM1’s stretched/unstretched if desired,
- computes thermodynamics from the familiar surface fields we are used to

This was done to allow controlled experiments with heterogeneous backgrounds to probe CM1 biases and investigate storm–environment feedbacks in controlled settings. If you have any questions (or beer offerings) feel free to reach out (tyler.pardun@noaa.gov)

---

## What’s in here

- base.F # Adds the isnd == 67 path to read and populate the CM1 grid with the 3D input_sounding
- cm1.F # Unchanged except accept adding xh and yh into the base.F subroutine
- cm1_init3d_nc.F # New NetCDF reader for 3D inputs
- cm1_interp3d.F # New horizontal and vertical resampling subroutines
- input_sounding.nc # This is the example input (python notebook on how to create it)
- Makefile # Modified build rules
- get_hrrr_for_cm1.ipynb #Notebook example to get a HRRR grid into the input_sounding.nc format

### New entry point
- `isnd = 67` → read 3-D analysis and build a heterogeneous base state.

---

## Required input file

A NetCDF file named `input_sounding.nc` with familiar "input_sounding" variables and their shapes:

- 3-D fields (stored as `(zh, yh, xh)` in the file):
  - `theta` [K], `qv` [g/kg], `u` [m/s], `v` [m/s]
- 2-D surface fields (stored as `(yh, xh)`):
  - `psfc` [hPa], `th0` [K], `qv0` [g/kg]
- Coordinates:
  - `xh(xh)`, `yh(yh)`, `zh(zh)` in **meters**, relative to some point origin

Example header:
```txt
dimensions:  zh = NZ ; yh = NY ; xh = NX
variables:
  float theta(zh,yh,xh) ;
  float qv(zh,yh,xh) ;
  float u(zh,yh,xh) ;
  float v(zh,yh,xh) ;
  float psfc(yh,xh) ;   // hPa
  float th0(yh,xh) ;    // K
  float qv0(yh,xh) ;    // g/kg
  double xh(xh) ; double yh(yh) ; double zh(zh) ; // meters
