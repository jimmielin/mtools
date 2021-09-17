# Plots

Tools for plotting. Generally follow the same data conventions and work out-of-the-box with data read by `Utils/Util_VarReader.m`.

Still preliminary.

* `Plot_DualCompare.m`: Plots a left/right/delta/ratio difference plot for a given species in a given level or time slice.
* `Plot_ZonalMean.m`: Plots zonal means.

## Data conventions
* Data is specified in `var(lon, lat, level, time)` (3-D data) or `var(lon, lat, time)` (2-D data) format.
* Level 1 is always assumed to be surface, with TOA facing "up" (indices increasing). **Util_VarReader.m** automatically rearranges data to fit this convention.
