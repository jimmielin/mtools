%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%            _____            ______        %
% _______ _____  /_______________  /_______ %
% __  __ `__ \  __/  __ \  __ \_  /__  ___/ %
% _  / / / / / /_ / /_/ / /_/ /  / _(__  )  %
% /_/ /_/ /_/\__/ \____/\____//_/  /____/   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       "mtools" Research Toolkit           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot_ZonalMean.m
%
% Plot zonal mean and diffs for two sets of data.
% (c) 2019-2021 Haipeng Lin <jimmie.lin@gmail.com>
%
% Version: 2021.09.15
% Started: 2021.09.15
% See changelog at end of script
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             PLOT PROPERTIES               %
%                                           %
%   <title> %s: species third%d fourth %d   %
%                                           %
%         LEFT | RIGHT | DELTA DIFF         %
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

% specify vertical limit for plot [hPa]
% tropopause is 70 hPa (~18km) ~ 400 hPa (~6km)
%
% usually, vert_top is set to 10 hPa
% surface, vert_bottom is set to 1000 hPa
vert_top = 10;
vert_bottom = 1000;

%%%%%%%%%%%%%%%%%%%%%% NO USER CONFIGURABLE CODE BELOW %%%%%%%%%%%%%%%%%%%%

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


% plot zonal mean given data, lat, lev, ASSUMING RECTILINEARITY
ax1 = subplot(1,3,1);

lats_plt = squeeze(lats_left(1,:));
levs_plt = squeeze(levs_left); % 0.999 sfc: 999 hPa conversion may be needed
data_plt_xyz = data_left(:,:,:,fourth_idx) .* convert_factor;

% collapse data by lats
data_plt_yz = squeeze(sum(data_plt_xyz)) / size(data_plt_xyz, 1); % sums over 1st dimension and avgd.
% this results in lats, levs which needs to be transposed
data_plt_zNy = flip(transpose(data_plt_yz), 1); % flip the vertical direction HACK so axes are TOA-up

if levs_plt(1) > levs_plt(2) % TOA is top level, flip the axis to make matlab happy
    levs_plt = flip(levs_plt);
end

contourf(lats_plt, levs_plt, data_plt_zNy, linspace(0, 10000, 10));
set(gca, 'YDir', 'reverse'); % ugly hack hack to get the axis labels right
set(gca, 'XDir', 'reverse'); % ugly hack hack to get the axis labels right
set(gca,'YScale', 'log');
colorbar;
colormap(ax1, flip(othercolor('BuDRd_12', 10)));
xlabel("Latitude [deg]");
ylabel("Pressure [hPa]");
title(left_name);

ylim([vert_top vert_bottom]); % 1 hPa is OK
    yticks([vert_top, flip(vert_bottom:-100:vert_top)]);
    xticks(linspace(-90, 90, 7));

set(gca, 'FontName', 'Arial');
set(gca, 'FontWeight', 'normal');
set(gca, 'FontSize', 14);
set(findall(gcf,'type','text'), 'FontName', 'Arial');
set(findall(gcf,'type','text'), 'FontWeight', 'normal');
set(findall(gcf,'type','text'), 'FontSize', 14);


%%%%%%%%%%
ax2 = subplot(1,3,2);

lats_plt = squeeze(lats_right(1,:));
levs_plt = squeeze(levs_right);
data_plt_xyz = data_right(:,:,:,fourth_idx) * convert_factor;

% collapse data by lats
data_plt_yz = squeeze(sum(data_plt_xyz)) / size(data_plt_xyz, 1); % sums over 1st dimension and avgd.
% this results in lats, levs which needs to be transposed
data_plt_zNy = flip(transpose(data_plt_yz), 1); % flip the vertical direction HACK so axes are TOA-up

if levs_plt(1) > levs_plt(2) % TOA is top level, flip the axis to make matlab happy
    levs_plt = flip(levs_plt);
end

contourf(lats_plt, levs_plt, data_plt_zNy, linspace(0, 10000, 10));
set(gca, 'YDir', 'reverse'); % ugly hack hack to get the axis labels right
set(gca, 'XDir', 'reverse'); % ugly hack hack to get the axis labels right
set(gca,'YScale', 'log');
colorbar;
colormap(ax2, flip(othercolor('BuDRd_12', 10)));
xlabel("Latitude [deg]");
ylabel("Pressure [hPa]");
title(right_name);


ylim([vert_top vert_bottom]); % 1 hPa is OK
    yticks([vert_top, flip(vert_bottom:-100:vert_top)]);
    xticks(linspace(-90, 90, 7));

% Set font properties...
set(gca, 'FontName', 'Arial');
set(gca, 'FontWeight', 'normal');
set(gca, 'FontSize', 14);
set(findall(gcf,'type','text'), 'FontName', 'Arial');
set(findall(gcf,'type','text'), 'FontWeight', 'normal');
set(findall(gcf,'type','text'), 'FontSize', 14);


%%%%%%%%%%%
% DO DELTA, IF POSSIBLE
if isequal(lats_left, lats_right) && isequal(levs_right, levs_left)
    ax2 = subplot(1,3,3);

    lats_plt = squeeze(lats_right(1,:));
    levs_plt = squeeze(levs_right);
    data_plt_xyz = (data_left(:,:,:,fourth_idx) - data_right(:,:,:,fourth_idx)) * convert_factor;

    % collapse data by lats
    data_plt_yz = squeeze(sum(data_plt_xyz)) / size(data_plt_xyz, 1); % sums over 1st dimension and avgd.
    % this results in lats, levs which needs to be transposed
    data_plt_zNy = flip(transpose(data_plt_yz), 1); % flip the vertical direction HACK so axes are TOA-up

    if levs_plt(1) > levs_plt(2) % TOA is top level, flip the axis to make matlab happy
        levs_plt = flip(levs_plt);
    end

    contourf(lats_plt, levs_plt, data_plt_zNy, linspace(-90, 90, 9));
    set(gca, 'YDir', 'reverse'); % ugly hack hack to get the axis labels right
    set(gca, 'XDir', 'reverse'); % ugly hack hack to get the axis labels right
    set(gca,'YScale', 'log');
    colorbar;
    colormap(ax2, flip(othercolor('RdBu11', 9)));
    xlabel("Latitude [deg]");
    ylabel("Pressure [hPa]");
    title("Delta [L-R]");
    
    caxis([-90 90]);

    ylim([vert_top vert_bottom]); % 1 hPa is OK
    yticks(flip(vert_bottom:-100:vert_top));
    xticks(linspace(-90, 90, 7));
    
    % Set font properties...
    set(gca, 'FontName', 'Arial');
    set(gca, 'FontWeight', 'normal');
    set(gca, 'FontSize', 14);
    set(findall(gcf,'type','text'), 'FontName', 'Arial');
    set(findall(gcf,'type','text'), 'FontWeight', 'normal');
    set(findall(gcf,'type','text'), 'FontSize', 14);
end

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

set(gcf, 'Renderer', 'painters', 'Position', [600 800 1800 580]);
