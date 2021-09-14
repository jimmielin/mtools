% read DAILY average data from HCO diagnostics outputs
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

filename = 'HEMCO_diagnostics.201605010000.nc';
prefixname = 'HEMCO_diagnostics.';
if strcmp(mode, 'daily')
    formatOut = 'yyyymmdd0000';
elseif strcmp(mode, 'hourly')
    formatOut = 'yyyymmddHHII';
end

suffixname = '.nc';

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
P0   = ncread(filename, 'P0');
% warning: this is only ref pressure as the surface pressure
% is not available in this diag?? FIXME hplin 3/18/21
levs = zeros(size(hyam, 1), 1);
for ll = 1:size(hyam,1)
    levs(ll) = P0 * hybm(ll) + hyam(ll);
    levs(ll) = levs(ll) / 100.0; % Pa to hPa
end
levs

Area_M2 = ncread(filename, 'AREA');

save("GC_HCO_coords.mat", "lons", "lats", "Area_M2");

%% actually read variables
% 3-D template:
spiename = 'EmisNO_Lightning';
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
EmisNO_Lightning = tmp(:,:,zdim:-1:(zdim-z_save+1),:);
disp(spiename);
save EmisNO_Lightning EmisNO_Lightning

% 2-D template
spiename = 'HcoLightningFlashRate_Total';
tmp = ncread(filename,spiename);
size(tmp)
xdim = size(tmp,1);
ydim = size(tmp,2);
tmp = zeros(xdim,ydim,size(datestd,1));
for tt = 1:size(datestd,1)
    filename = [prefixname,date_str(tt,:),suffixname];
    disp(filename);
    tmp(:,:,tt) = ncread(filename,spiename);
end
HcoLightningFlashRate_Total = tmp(:,:,:);
disp(spiename);
save HcoLightningFlashRate_Total HcoLightningFlashRate_Total

%% collapse EmisNO_Lightning into column
EmisNO_Lightning_Col = zeros(xdim,ydim,size(datestd,1));
for i = 1:xdim
    for j = 1:ydim
        for tt = 1:size(datestd,1)
            EmisNO_Lightning_Col(i,j,tt) = sum(EmisNO_Lightning(i,j,:,tt));
        end
    end
end
save EmisNO_Lightning_Col EmisNO_Lightning_Col
