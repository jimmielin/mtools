% read DAILY average data from CESM outputs
% hplin, 3/6/21
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

filename = 'case_cesm-gc_feb21.cam.h1.2016-05-01-00000.nc';
prefixname = 'case_cesm-gc_feb21.cam.h1.';
suffixname = '.nc';

% z-levels to save? (for 3-D data)
% SAVED FROM BOTTOM-UP, even though TOA is lev #1.
% after reading, data will be inverted (TOA is lev #max)
z_save = 56;

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
formatOut = 'yyyy-mm-dd-00000';
date_str = datestr(datestd,formatOut);

%% read the coords
lons = ncread(filename, 'lon');
lats = ncread(filename, 'lat');
levs = ncread(filename, 'lev');
levs = flip(levs); % this is in Pa for CESM format

if size(levs, 1) < z_save
    fprintf("WARNING: z_save > max levs. Setting as max available.\n")
    z_save = size(levs, 1)
end

save("CESM_coords.mat", "lons", "lats", "levs");

%% actually read variables

% 3-D template:
spiename = 'O3';
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
O3 = tmp(:,:,zdim:-1:(zdim-z_save+1),:);
disp(spiename);
save O3 O3
