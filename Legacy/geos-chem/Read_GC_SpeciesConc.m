% read DAILY average data from GEOS-Chem SpeciesConc outputs
% hplin, 3/16/21
% based off original WRF-GC script by xfeng

mode           = 'daily';

rundays        = 7;
startyr        = 2016;
endyr          = 2016;
startmon       = 5;
endmon         = 5;
startdate      = 1; 
enddate        = 30;
num_days       = datenum(endyr,endmon,enddate) - datenum(startyr,startmon,startdate) + 1;

filename = 'GEOSChem.SpeciesConc.20160501_0000z.nc4';
prefixname = 'GEOSChem.SpeciesConc.';
if strcmp(mode, 'daily')
    formatOut = 'yyyymmdd_0000';
elseif strcmp(mode, 'hourly')
    formatOut = 'yyyymmdd_HHII';
end

suffixname = 'z.nc4';

% z-levels to save? (for 3-D data)
% SAVED FROM BOTTOM-UP, even though TOA is lev #1.
% after reading, data will be inverted (TOA is lev #max)
z_save = 72;

% non-editable code below

if strcmp(mode, 'daily')
    datestd        = zeros(num_days, 6);
elseif strcmp(mode, 'hourly')
    datestd        = zeros(24 * num_days, 6);
    % datestd dims: Y,M,D,H,I,S
end
        
i = 0;
for mm = datenum(startyr,startmon,startdate):datenum(endyr,endmon,enddate)
    if strcmp(mode, 'hourly')
        for hh = 0:23
            i = i + 1;
            datestd(i,:) = datevec(mm);
            datestd(i,4) = hh;
        end
    elseif strcmp(mode, 'daily')
        i = i + 1;
        datestd(i,:) = datevec(mm);
        datestd(i,4) = 0; % always 0-hour. change for slicing
    end
end

% formatOut = 'yyyy-mm-dd_HH:MM:SS';

% assumes 1 slice per file, 1 file per day
% (can also be 24 slices if using 'hourly')
date_str = datestr(datestd,formatOut);

%% read lat, lon, lev for storage
lons = ncread(filename, 'lon');
lats = ncread(filename, 'lat');

% process into pressure-level levs
hyam = ncread(filename, 'hyam');
hybm = ncread(filename, 'hybm');
levs = ncread(filename, 'lev');
levs

save("GC_GC_coords.mat", "lons", "lats", "levs");

%% actually read variables
% 3-D template:
spiename = 'SpeciesConc_O3';
tmp = ncread(filename,spiename);

% in CESM, data is time,lev,lat,lon ? check
% 288x192x56 lon,lat,lev so just read like WRF
size(tmp)
xdim = size(tmp,1);
ydim = size(tmp,2);
zdim = size(tmp,3);
tmp = zeros(xdim,ydim,zdim,size(datestd,1));
for tt = 1:size(datestd,1)
    filename = [prefixname,date_str(tt,:),suffixname];
    disp(filename);
    tmp(:,:,:,tt) = ncread(filename,spiename);
end
GC_O3 = tmp(:,:,1:z_save,:);
disp(spiename);
save GC_O3 GC_O3

% 2-D template
