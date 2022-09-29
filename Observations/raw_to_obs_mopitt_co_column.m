% get data: https://search.earthdata.nasa.gov/search/granules?p=C1575974140-LARC&pg[0][v]=f&pg[0][qt]=2016-01-01T00%3A00%3A00.000Z%2C2017-01-01T23%3A59%3A59.999Z&pg[0][gsk]=-start_date&g=G1576022220-LARC&tl=1663289553!3!!
% MOPITT CO gridded monthly means (Near and Thermal Infrared Radiances)V008

data_source_fstring = "raw_mopitt_co/MOP03JM-%04d%02d-L3V95.6.3.he5";
data_sink_fstring = "obs_mopitt_co/%04d%02d_monthly_mean.mat";

year = 2016;
month = 8;

data_source = sprintf(data_source_fstring, year, month);
data_sink = sprintf(data_sink_fstring, year, month);

mopitt_co_column_day = h5read(data_source, '/HDFEOS/GRIDS/MOP03/Data Fields/RetrievedCOTotalColumnDay');
mopitt_co_column_night = h5read(data_source, '/HDFEOS/GRIDS/MOP03/Data Fields/RetrievedCOTotalColumnNight');

mopitt_co_sfc_ppb_day = h5read(data_source, '/HDFEOS/GRIDS/MOP03/Data Fields/RetrievedCOSurfaceMixingRatioDay');
mopitt_co_sfc_ppb_night = h5read(data_source, '/HDFEOS/GRIDS/MOP03/Data Fields/RetrievedCOSurfaceMixingRatioNight');

mopitt_co_column_day(mopitt_co_column_day < 0) = NaN;
mopitt_co_column_night(mopitt_co_column_night < 0) = NaN;
mopitt_co_sfc_ppb_day(mopitt_co_sfc_ppb_day < 0) = NaN;
mopitt_co_sfc_ppb_night(mopitt_co_sfc_ppb_night < 0) = NaN;

mopitt_lons = h5read(data_source, '/HDFEOS/GRIDS/MOP03/XDim');
mopitt_lats = h5read(data_source, '/HDFEOS/GRIDS/MOP03/YDim');

save(data_sink, "mopitt_co_column_day", "mopitt_co_column_night", "mopitt_co_sfc_ppb_day", "mopitt_co_sfc_ppb_night", "mopitt_lons", "mopitt_lats");