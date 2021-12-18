# Metrics

Tools for computing data diagnostics. Generally follow the same data conventions and work out-of-the-box with data read by `Utils/Util_VarReader.m`.

Still preliminary.

* `Metrics_GlobalMass.m` computes two-model global mass of given species and percentage diffs. It also plots column ozone if `O3` variable is requested.

## Data conventions
* Data is specified in `var(lon, lat, level, time)` (3-D data) or `var(lon, lat, time)` (2-D data) format.
* Level 1 is always assumed to be surface, with TOA facing "up" (indices increasing). **Util_VarReader.m** automatically rearranges data to fit this convention.

## Tools which do not confirm to specifications in `mtools` for simplicity
* `Metrics_Diff_RMSE.m` computes **two-run** (**not** two-model) RRMS and RMSEs of given species, vertical level, and whole species array. It can also create RRMS histograms, by-level plots, vertically integrated plots by species, by-species totals, and totals. It can do geographical subsetting and be forced from an external time step to save histogram peaks of RRMS by-species. The data is read directly from raw netCDF files for convenience.
