%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%            _____            ______        %
% _______ _____  /_______________  /_______ %
% __  __ `__ \  __/  __ \  __ \_  /__  ___/ %
% _  / / / / / /_ / /_/ / /_/ /  / _(__  )  %
% /_/ /_/ /_/\__/ \____/\____//_/  /____/   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       "mtools" Research Toolkit           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot_DualCompare.m
%
% Plot arbitrary comparison plots between *two* sets of data.
% with the option to do a delta and ratio diff.
% (c) 2019-2021 Haipeng Lin <jimmie.lin@gmail.com>
%
% Version: 2021.09.14
% Started: 2020.06.08
% See changelog at end of script
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             PLOT PROPERTIES               %
%                                           %
%   <title> %s: species third%d fourth %d   %
%                                           %
%    L  E  F  T     |   R  I  G  H  T       %
%                                           %
%     DELTA DIFF    |     RATIO DIFF        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% specify data coordinates.
%
% if data is not on the same grid, they will be regridded.
% provide coordinates in following form:

coords_CESM = load('coords_CESM.mat');
lons_left   = coords_CESM.lons;
lats_left   = coords_CESM.lats;
levs_left   = coords_CESM.levs;

lons_right  = coords_CESM.lons;
lats_right  = coords_CESM.lats;
levs_right  = coords_CESM.levs;

% special data coordinate handling.
% some data is fancy, and does not fit the m_map specs.
% so some processing has to be done, such as circshifting data, or expand
% coordinates from 1-D arrays to 2-D.
%
% specify the handling below with the following magic values:
% =>   cesm-circshift: coordinates come from CESM and are shifted 180deg.
%                      both DATA and COORDS will be shifted to fix.
% =>   gc:             coordinates are in 1-D and need to be changed to 2D
% everything else is ignored.
handling_left  = "cesm-circshift";
handling_right = "cesm-circshift";

% specify data on left and right
%
% data needs to be 4-D, with dimensions (lon, lat, third_idx, fourth_idx)
% where third_idx is /usually/ the level, and
%       fourth_idx is /usually/ a date slice.
%
% 
file = load('climo_CAMchemYearlyMonthlyAvg_2016.mat');
data_left = file.yearly_O3;

file = load('climo_CESM2GCYearlyMonthlyAvg_2016.mat');
data_right = file.O3;

% specify plot properties
%
% title:       format string where %s will insert the species name (var_name)
%              and %d, %d will insert:
%                 the third index of the data & the fourth index.
% var_name:    the "nice" name of the variable, including units
% left_name:   source of left data
% right_name:  source of right data
title_nice     = "2016 monthly average - %s - level %d - date slice %d";
var_nice_name  = "O3 [ppb]";

left_name      = "CAM-chem";
right_name     = "CESM2-GC";

% specify data indexing properties
%
% third_idx: /usually/ level. for CESM 56-levels, 23 is ~500 hPa.
%            magic number -1 will make it "column sum" (summed over index)
% fourth_idx: usually date slice
%
% WARNING: IF DATA IS GIVEN IN 3-DIMENSIONAL ARRAY, USUALLY FOR 2-D DATA
% IT WILL BE INTERPRETED AS (lon, lat, fourth_idx)
% BECAUSE IT USUALLY MEANS IT IS 2-D DATA, AND FOURTH_IDX IS DATE SLICE BY
% CONVENTION.
third_idx =  1; % level index
fourth_idx = 4; % date index

% specify data scaling properties
%
% usually for ozone, 1e9 to convert to ppb.
% otherwise, 1
convert_factor = 1e9;

% specify data range properties
%
% we do not make educated guesses here, so you have to give your own.
% there are some representative sizes below
%
% maxDynRange: max dynamic range for value
% maxDynRangeCmp: max dynamic range for comparison (delta) plot
% distRatioCmp: dynamic range half-size for ratio plot
%               e.g. 0.25 = ratio plot from 0.75 ~ 1.25
%
%      level |    surface    |      500 hPa     |     column sum    |
%   /species |------------------------------------------------------|
% -----------|                                                      |
%      O3    | 60 ppb, d. 18 |  90 ppb, d. 20   |   2.4e5, d. 3e4   |
%     NO2    | 10 ppb, d.4.5 |
%   J-O3O1D  |  0.000012,d/2 |
%   J-O3O3P  |  0.0006  ,d/2 |
%   J-NO2    |  0.006   ,d/2 |
%   J-CH2O   |  0.00012,     |
%            |  d. 000045    |
% -------------------------------------------------------------------

maxDynRange = 60;
maxDynRangeCmp = 18;
distRatioCmp = 0.5;

% maxDynRange = 12;
% maxDynRangeCmp = 4.5;
% distRatioCmp = 0.5;




%%%%%%%%%%%%%%%%%%%%%% NO USER CONFIGURABLE CODE BELOW %%%%%%%%%%%%%%%%%%%%

% handling of coordinates ...
% cesm-circshift
if strcmp(handling_left, 'cesm-circshift')
    % for CESM processing - if lons are at 0.0 to 360.0, need to transpose data
    % into -180.0, 180.0. check
    if max(lons_left) > 180.0
       data_left = circshift(data_left, ...
                     size(lons_left, 1)/2, ...
                     1); % data needs to be in (lat, lon) 

       lons_left = lons_left - 180.0;

       fprintf("cesm-circshift: shifted data_left, lons_left\n");
    end
end

if strcmp(handling_right, 'cesm-circshift')
    % for CESM processing - if lons are at 0.0 to 360.0, need to transpose data
    % into -180.0, 180.0. check
    if max(lons_right) > 180.0
       data_right = circshift(data_right, ...
                     size(lons_right, 1)/2, ...
                     1); % data needs to be in (lat, lon) 

       lons_right = lons_right - 180.0;

       fprintf("cesm-circshift: shifted data_right, lons_right\n");
    end
end

% auto-fix 2-D coordinates
if size(lons_left, 2) == 1
    lons_new = zeros(size(lons_left, 1), size(lats_left, 1));
    lats_new = zeros(size(lons_left, 1), size(lats_left, 1));

    for j = 1:size(lats_left, 1)
        lons_new(:,j) = lons_left(:,1);
    end

    for i = 1:size(lons_left, 1)
        lats_new(i,:) = lats_left(:,1);
    end

    lons_left = lons_new;
    lats_left = lats_new;
end

if size(lons_right, 2) == 1
    lons_new = zeros(size(lons_right, 1), size(lats_right, 1));
    lats_new = zeros(size(lons_right, 1), size(lats_right, 1));

    for j = 1:size(lats_right, 1)
        lons_new(:,j) = lons_right(:,1);
    end

    for i = 1:size(lons_right, 1)
        lats_new(i,:) = lats_right(:,1);
    end

    lons_right = lons_new;
    lats_right = lats_new;
end

% HANDLE INDEX SETTING
if length(size(data_left)) == 4
    if third_idx > 0
        data_left  = data_left (:,:,third_idx,fourth_idx) * convert_factor;
        data_right = data_right(:,:,third_idx,fourth_idx) * convert_factor;
    else
        data_left  = sum(data_left (:,:,:,fourth_idx), 3) * convert_factor;
        data_right = sum(data_right(:,:,:,fourth_idx), 3) * convert_factor;
    end
elseif length(size(data_left)) == 3
    fprintf("Note: 2-D data detected; data is interpreted as lon,lat,fourth_idx. buyer beware.\n")
    
    data_left  = data_left (:,:,fourth_idx) * convert_factor;
    data_right = data_right(:,:,fourth_idx) * convert_factor;
else
    error("size of data_left is neither 3-D or 4-D. what is going on?");
end

% CHECK IF REGRIDDING IS NECESSARY.
if isequal(lats_left, lats_right) && isequal(lons_left, lons_right)
    data_diff = data_left - data_right;
    data_ratio = data_left ./ data_right;
else
    % todo: we force regridding to right grid
    data_l2r = griddata(lons_left, lats_left, data_left, ...
                        lons_right, lats_right);
    data_diff = data_l2r - data_right;
    data_ratio = data_l2r ./ data_right;
    
    % for comparison plots, right-coordinate is always used for now.
    % (hplin, 9/14/21)
end

% MAKE subplot_tightS
% LEFT
axg = subplot_tight(2,2,1);

hold on;
m_proj('miller', ...
       'long',[-180.0 180.0], ...
       'lat', [-87.5 87.5]);

m_pcolor(lons_left, lats_left, data_left);
colorbar;
caxis([0, maxDynRange]);

m_plotbndry('/usr/local/MATLAB/R2021a/toolbox/boundary/world', ...
    'color', 'k', 'linewidth', 1.5) % range lon is -180.0 180.0 so correct
%m_plotbndry('C:\Program Files\MATLAB\R2021a\toolbox\boundary\shengjie', ...
%    'color', 'k', 'linewidth', 1.5) % range lon is -180.0 180.0 so correct
m_grid;
colormap(axg, flip(othercolor('Spectral6', 12)));

t = title(sprintf(left_name));

% RIGHT
axg = subplot_tight(2,2,2);

hold on;
m_proj('miller', ...
       'long',[-180.0 180.0], ...
       'lat', [-87.5 87.5]);

m_pcolor(lons_right, lats_right, data_right);
colorbar;
caxis([0, maxDynRange]);

m_plotbndry('/usr/local/MATLAB/R2021a/toolbox/boundary/world', ...
    'color', 'k', 'linewidth', 1.5) % range lon is -180.0 180.0 so correct
%m_plotbndry('C:\Program Files\MATLAB\R2021a\toolbox\boundary\shengjie', ...
%    'color', 'k', 'linewidth', 1.5) % range lon is -180.0 180.0 so correct
m_grid;
colormap(axg, flip(othercolor('Spectral6', 12)));

t = title(sprintf(right_name));

% DELTAS
ax3 = subplot_tight(2,2,3);
m_proj('miller', ...
       'long',[-180.0 180.0], ...
       'lat', [-87.5 87.5]);

m_pcolor(lons_right, lats_right, data_diff);
hold on;
colorbar;

caxis([-maxDynRangeCmp, maxDynRangeCmp]);

m_plotbndry('/usr/local/MATLAB/R2021a/toolbox/boundary/world', ...
    'color', 'k', 'linewidth', 1.5)
m_grid;
colormap(ax3, flip(othercolor('RdBu11', 9)));

t = title(sprintf("Diff (L-R)"));

% RATIO
ax4 = subplot_tight(2,2,4);
m_proj('miller', ...
       'long',[-180.0 180.0], ...
       'lat', [-87.5 87.5]);

m_pcolor(lons_right, lats_right, data_ratio);
hold on;
colorbar;

caxis([1 - distRatioCmp, 1 + distRatioCmp]);

m_plotbndry('/usr/local/MATLAB/R2021a/toolbox/boundary/world', ...
    'color', 'k', 'linewidth', 1.5)
m_grid;
colormap(ax4, flip(othercolor('RdBu11', 9)));

t = title(sprintf("Ratio (L/R)"));

sgtitle(sprintf(title_nice, var_nice_name, third_idx, fourth_idx));


% Set font properties...
set(gca, 'FontName', 'Arial');
set(gca, 'FontWeight', 'normal');
set(gca, 'FontSize', 14);
set(findall(gcf,'type','text'), 'FontName', 'Arial');
set(findall(gcf,'type','text'), 'FontWeight', 'normal');
set(findall(gcf,'type','text'), 'FontSize', 14);

% Set printing properties
%set(gcf, 'PaperUnits', 'inches', 'PaperPosition', [0 0 12 10]);
%set(gcf, 'Renderer', 'painters', 'Position', [90 90 1800 500])

set(gcf, 'Renderer', 'painters', 'Position', [90 0 1250 970]);
