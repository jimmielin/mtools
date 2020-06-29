% Read data from GEOSChem SpeciesConc netCDF Diagnostics
% and save to mat file for later use
% Haipeng Lin, 2020/06/08

fname_prefix = 'GEOSChem.SpeciesConc.';
fname_suffix = 'z.nc4';

s_year  = 2016;
s_month = 7;
s_day   = 1;
e_year  = 2016;
e_month = 7;
e_day   = 4; % including this day

delta_h = 1;

%% Create GC Dates
num_days = datenum(e_year,e_month,e_day) - ...
           datenum(s_year,s_month,s_day) + 1;
datestd  = zeros(24/delta_h * num_days, 6); % 2nd dim is YMDHIS

i = 0;
for mm = datenum(s_year,s_month,s_day):datenum(e_year,e_month,e_day)
    for hh = 0:delta_h:(24-delta_h)
        i = i + 1;
        datestd(i,:) = datevec(mm);
        datestd(i,4) = hh;
    end
end

formatGC = 'yyyymmdd_HHMM';
date_strGC = datestr(datestd(1:end-3,:), formatGC);

numtimes = datenum(datestd(1:end-3,:)) - 1;

%% Read the lat-lon parameters
fname = [fname_prefix, date_strGC(1,:), fname_suffix];

lats = ncread(fname, 'lat');
lons = ncread(fname, 'lon');
levs = ncread(fname, 'lev');
ftimes = ncread(fname, 'time');
times = length(date_strGC); % fixme: assume 1 time per file

%% Only plotting and saving surface data for now (lev = 1)
tmp2D_NH3 = zeros(length(lons), length(lats), times);
tmp2D_NO  = zeros(length(lons), length(lats), times);
tmp2D_NO2 = zeros(length(lons), length(lats), times);
tmp2D_O3 = zeros(length(lons), length(lats), times);
tmp2D_SO2 = zeros(length(lons), length(lats), times);
tmp4D_read = zeros(length(lons), length(lats), length(levs), length(ftimes));

for t = 1:length(date_strGC)
    fname = [fname_prefix, date_strGC(t,:), fname_suffix];
    disp(fname);
    
    if not(isfile(fname))
        disp(strcat("Could not find file: ", fname));
        break
    end
    
    tmp4D_read = ncread(fname, 'SpeciesConc_NH3');
    tmp2D_NH3(:,:,t) = tmp4D_read(:,:,1,1); % lev, time collapsed
    
    tmp4D_read = ncread(fname, 'SpeciesConc_NO');
    tmp2D_NO(:,:,t) = tmp4D_read(:,:,1,1); % lev, time collapsed
    
    tmp4D_read = ncread(fname, 'SpeciesConc_NO2');
    tmp2D_NO2(:,:,t) = tmp4D_read(:,:,1,1); % lev, time collapsed
    
    tmp4D_read = ncread(fname, 'SpeciesConc_O3');
    tmp2D_O3(:,:,t) = tmp4D_read(:,:,1,1); % lev, time collapsed
    
    tmp4D_read = ncread(fname, 'SpeciesConc_SO2');
    tmp2D_SO2(:,:,t) = tmp4D_read(:,:,1,1); % lev, time collapsed
end

save('tmp_speciesconc_2.mat', 'lats', 'lons', 'numtimes', ...
     'tmp2D_NH3', 'tmp2D_NO', 'tmp2D_NO2', 'tmp2D_O3', 'tmp2D_SO2');