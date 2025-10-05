# CM1 Heterogeneous Initialization (isnd=67)

Modify CM1 to initialize from a **3-D heterogeneous analysis** (e.g., HRRR analysis). This adds a new initialization path `isnd=67` that:
- reads a NetCDF file with 3-D fields on a source grid,
- horizontally resamples each source level to the CM1 horizontal grid,
- vertically and horizontally interpolates to CM1’s stretched/unstretched,
- computes thermodynamics hydrostatically from the provided surface fields,
- populates CM1 base arrays (`th0, qv0, u0, v0, pi0, prs0, t0, rh0`) on the mass grid.

**Goal:** allow controlled experiments with heterogeneous backgrounds to probe CM1 biases and storm–environment feedbacks.

---

## What’s in here

base.F # adds ELSEIF (isnd == 67) path
cm1.F # main driver (unchanged except accept adding xh and yh into the base subroutine)
cm1_init3d_nc.F # read_init3d_nc(): NetCDF reader for 3D inputs
cm1_interp3d.F # bilinear_resample2d() and vert_linear_to_zh()
input_sounding.nc # example input (layout below)
Makefile # Modified build rules
get_hrrr_for_cm1.ipynb #Notebook example to get a HRRR grid into the input_sounding.nc format

### New entry point
- `isnd = 67` → read 3-D analysis and build a heterogeneous base state.
- Leaves `isnd = 7` behavior intact (homogeneous tiling of 1-D profile).

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
