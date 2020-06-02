opengl hardware;
colormap(flip(othercolor('RdYlBu10', 8)));

% Load netCDF data (UTC time)
fname = 'case_cesm_gc.cam.h1.2000-01-02-00000.nc';
hist_ncid = netcdf.open(fname, 'NC_NOWRITE');
hist_varid_lons = netcdf.inqVarID(hist_ncid, 'lon');
hist_varid_lats = netcdf.inqVarID(hist_ncid, 'lat');
hist_varid_levs = netcdf.inqVarID(hist_ncid, 'lev');
hist_varid_times = netcdf.inqVarID(hist_ncid, 'time');
hist_varid_p0 = netcdf.inqVarID(hist_ncid, 'P0');
hist_varid_hyam = netcdf.inqVarID(hist_ncid, 'hyam');
hist_varid_hybm = netcdf.inqVarID(hist_ncid, 'hybm');

hist_lons = netcdf.getVar(hist_ncid, hist_varid_lons);
hist_lats = netcdf.getVar(hist_ncid, hist_varid_lats);
hist_levs = netcdf.getVar(hist_ncid, hist_varid_levs);
hist_p0   = netcdf.getVar(hist_ncid, hist_varid_p0);
hist_hyams = netcdf.getVar(hist_ncid, hist_varid_hyam);
hist_hybms = netcdf.getVar(hist_ncid, hist_varid_hybm);
hist_times = netcdf.getVar(hist_ncid, hist_varid_times);

varname = 'HCO_NO2';
hist_varid_hcovar = netcdf.inqVarID(hist_ncid, varname);
hist_hcovar = netcdf.getVar(hist_ncid, hist_varid_hcovar);
% hist_hcovar_meta = ncinfo(fname, varname); % not working well atm


% if hist_times is only 1 dimension, then hist_hcovar will be lonxlatxlevx1
% collapsed to 3-d only

% for reusable code
lons = double(hist_lons(:,1)) - 180.0; % orig. 0 to 360
lats = double(hist_lats(:,1));
lev  = size(hist_levs, 1);
plev = hist_hyams(lev) + hist_hybms(lev) * hist_p0/100; % hPa
data = circshift(transpose(hist_hcovar(:,:,lev,1)), ...
                 size(hist_lons, 1)/2, ...
                 2); % data needs to be in (lat, lon)

% fill in the periodic boundary (1peridim)
lons = [lons; 180.0];    % vertcat
data = [data data(:,1)]; % horzcat

% Set display properties
set(gcf, 'Renderer', 'painters', 'Position', [90 90 1000 600])

% [min(min(lons)) max(max(lons))]
% [min(min(lats)) max(max(lats))]
m_proj('miller', ...
       'long',[-180.0 180.0], ...
       'lat', [-89.95 89.95]);

m_pcolor(lons, lats, data); % data is in (lat, lon) whereas lons x lats input...
hold on;
colorbar;

m_plotbndry('C:\Program Files\MATLAB\R2020a\toolbox\boundary\world', ...
    'color', 'k', 'linewidth', 1.5) % range lon is -180.0 180.0 so correct
m_grid;

t = title(strcat("HEMCO Emissions Flux ", strrep(varname, 'HCO_', ''), ...
                 " Level ", int2str(lev), ' (', num2str(plev), ' hPa) [kg m^{2} s^{-1}]'));

% Set font properties...
set(gca, 'FontName', 'Arial');
set(gca, 'FontWeight', 'normal');
set(gca, 'FontSize', 14);
set(findall(gcf,'type','text'), 'FontName', 'Arial');
set(findall(gcf,'type','text'), 'FontWeight', 'normal');
set(findall(gcf,'type','text'), 'FontSize', 14);

% Set printing properties
set(gcf, 'PaperUnits', 'inches', 'PaperPosition', [0 0 12 10]);