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
% (c) 2019-2022 Haipeng Lin <jimmie.lin@gmail.com>
%
% Version: 2022.09.28
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

coords_f19 = load('mdl_coords_CESM_f19_f19_mg17.mat');
coords_f09 = load('mdl_coords_CESM_f09_f09_mg17.mat');
coords_left  = coords_f19;
coords_right = coords_f19;

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
fname_left  = "mdl_cesm2.1gc13.4.1_2016_monthly.mat";
%fname_right = "mdl_cesm2.2camchem_2016_monthly.mat";
fname_right = "mdl_cesm2.1gc13.1.2_2016_monthly.mat";

var_name = "O3";

file_left = load(fname_left, var_name);
data_left = file_left.(var_name);
file_right = load(fname_right, var_name);
data_right = file_right.(var_name);

%% non-data intensive section
% specify data scaling properties
%
% usually for ozone, 1e9 to convert to ppb.
% otherwise, 1
convert_factor = 1e9;
if convert_factor == 1
    var_unit = "VMR";
elseif convert_factor == 1e9
    var_unit = "ppbv";
elseif convert_factor == 1e12
    var_unit = "pptv";
else
    var_unit = sprintf("%e VMR", convert_factor);
end

% specify plot properties
%
% title:       format string where %s will insert the species name (var_name)
%              and %d, %d will insert:
%                 the third index of the data & the fourth index.
% var_name:    the "nice" name of the variable, including units
% left_name:   source of left data
% right_name:  source of right data
title_nice     = "2016 monthly average - %s - level %d - date slice %d";
var_nice_name  = sprintf("%s [%s]", var_name, var_unit);

left_name      = "CESM2.1-GC13.4.1";
%right_name     = "CESM2.2-CAMChem";
right_name = "CESM2.1-GC13.1.2";

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
fourth_idx = 7; % date index

% specify vertical limit for plot [hPa]
% tropopause is 70 hPa (~18km) ~ 400 hPa (~6km)
%
% usually, vert_top is set to 10 hPa
% surface, vert_bottom is set to 1000 hPa
vert_top = 100;
vert_bottom = 1000;

% max dynamic range
% ozone: 10, 100, 30 [ppb]
% OH: 0, 0.4, 0.2 [ppt]
minDynRange = 10;
maxDynRange = 100;
maxDynRangeCmp = 45;
distRatioCmp = 1;

%%%%%%%%%%%%%%%%%%%%%% NO USER CONFIGURABLE CODE BELOW %%%%%%%%%%%%%%%%%%%%
% ticks
pressure_ticks = [vert_top, flip(vert_bottom:-100:vert_top)];
if pressure_ticks(1) == pressure_ticks(2)
    pressure_ticks = pressure_ticks(2:end);
end

% coords
lons_left   = coords_left.lons;
lats_left   = coords_left.lats;
levs_left   = coords_left.levs;

lons_right  = coords_right.lons;
lats_right  = coords_right.lats;
levs_right  = coords_right.levs;

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

% Compute dynamic range for yTicks calculation.
% Special handling if the maximum value is < 1, as there will be scaling
if maxDynRange > 1
    yTicksAbs = minDynRange:ceil(maxDynRange/10):maxDynRange;
else
    % try scale_factor that works
    yTicks_scaleFactor = 10;
    dynRangeTest = maxDynRange * yTicks_scaleFactor;
    while yTicks_scaleFactor < 1e16 && dynRangeTest < 1
        yTicks_scaleFactor = yTicks_scaleFactor * 10;
        dynRangeTest = maxDynRange * yTicks_scaleFactor;
    end
    
    if dynRangeTest < 1
        error("Are values too small? We cannot find an appropriate scale.");
    else
        % settle on a scale
        yTicksAbs = (minDynRange*yTicks_scaleFactor*10):ceil(dynRangeTest):dynRangeTest*10;
        yTicksAbs = yTicksAbs / yTicks_scaleFactor / 10;
    end
end

% plot zonal mean given data, lat, lev, ASSUMING RECTILINEARITY
ax1 = subplot(1,4,1);

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

contourf(lats_plt, levs_plt, data_plt_zNy, linspace(minDynRange, maxDynRange, 10));
set(gca, 'YDir', 'reverse'); % ugly hack hack to get the axis labels right
set(gca, 'XDir', 'reverse'); % ugly hack hack to get the axis labels right
set(gca,'YScale', 'log');
caxis([minDynRange maxDynRange]);
colorbar;
colormap(ax1, (othercolor('YlGnBu9', 9)));
xlabel("Latitude [deg]");
ylabel("Pressure [hPa]");
title(left_name);

ylim([vert_top vert_bottom]); % 1 hPa is OK
yticks(pressure_ticks);
xticks(linspace(-90, 90, 7));

set(gca, 'FontName', 'Arial');
set(gca, 'FontSize', 16);
set(findall(gcf,'type','text'), 'FontName', 'Arial');
set(findall(gcf,'type','text'), 'FontSize', 16);


%%%%%%%%%%
ax2 = subplot(1,4,2);

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

contourf(lats_plt, levs_plt, data_plt_zNy, linspace(minDynRange, maxDynRange, 10));
set(gca, 'YDir', 'reverse'); % ugly hack hack to get the axis labels right
set(gca, 'XDir', 'reverse'); % ugly hack hack to get the axis labels right
set(gca,'YScale', 'log');
caxis([minDynRange maxDynRange]);
colorbar;
colormap(ax2, (othercolor('YlGnBu9', 9)));
xlabel("Latitude [deg]");
%ylabel("Pressure [hPa]");
title(right_name);


ylim([vert_top vert_bottom]); % 1 hPa is OK
yticks(pressure_ticks);
xticks(linspace(-90, 90, 7));

% Set font properties...
set(gca, 'FontName', 'Arial');
set(gca, 'FontSize', 16);
set(findall(gcf,'type','text'), 'FontName', 'Arial');
set(findall(gcf,'type','text'), 'FontSize', 16);


%%%%%%%%%%%
% DO DELTA, IF POSSIBLE
if isequal(lats_left, lats_right) && isequal(levs_right, levs_left)
    ax2 = subplot(1,4,3);

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

    contourf(lats_plt, levs_plt, data_plt_zNy, linspace(-maxDynRangeCmp, maxDynRangeCmp, 9));
    set(gca, 'YDir', 'reverse'); % ugly hack hack to get the axis labels right
    set(gca, 'XDir', 'reverse'); % ugly hack hack to get the axis labels right
    set(gca,'YScale', 'log');
    cbrd = colorbar;
    colormap(ax2, flip(othercolor('RdBu11', 9)));
    xlabel("Latitude [deg]");
    %ylabel("Pressure [hPa]");
    title("Delta [L-R]");
    
    set(cbrd, 'YTick', round(unique([-maxDynRangeCmp:maxDynRangeCmp/4.5:0 0 maxDynRangeCmp/9:maxDynRangeCmp/4.5:maxDynRangeCmp]), 2));
    caxis([-maxDynRangeCmp maxDynRangeCmp]);

    ylim([vert_top vert_bottom]); % 1 hPa is OK
    yticks(flip(vert_bottom:-100:vert_top));
    xticks(linspace(-90, 90, 7));
    
    % Set font properties...
    set(gca, 'FontName', 'Arial');
    set(gca, 'FontSize', 16);
    set(findall(gcf,'type','text'), 'FontName', 'Arial');
    set(findall(gcf,'type','text'), 'FontSize', 16);
    
    
    ax2 = subplot(1,4,4);

    lats_plt = squeeze(lats_right(1,:));
    levs_plt = squeeze(levs_right);
    data_plt_xyz = (data_left(:,:,:,fourth_idx) ./ data_right(:,:,:,fourth_idx));

    % collapse data by lats
    data_plt_yz = squeeze(sum(data_plt_xyz)) / size(data_plt_xyz, 1); % sums over 1st dimension and avgd.
    % this results in lats, levs which needs to be transposed
    data_plt_zNy = flip(transpose(data_plt_yz), 1); % flip the vertical direction HACK so axes are TOA-up

    if levs_plt(1) > levs_plt(2) % TOA is top level, flip the axis to make matlab happy
        levs_plt = flip(levs_plt);
    end

    contourf(lats_plt, levs_plt, data_plt_zNy, linspace(1-distRatioCmp, 1+distRatioCmp, 10));
    set(gca, 'YDir', 'reverse'); % ugly hack hack to get the axis labels right
    set(gca, 'XDir', 'reverse'); % ugly hack hack to get the axis labels right
    set(gca,'YScale', 'log');
    cbrr = colorbar;
    colormap(ax2, flip(othercolor('RdBu11', 9)));
    xlabel("Latitude [deg]");
    %ylabel("Pressure [hPa]");
    title("Ratio [L/R]");
    
    set(cbrr, 'YTick', round(unique([1-distRatioCmp:distRatioCmp/4.5:1 1 (1+distRatioCmp/9):distRatioCmp/4.5:1+distRatioCmp]), 2));
    
    caxis([1-distRatioCmp 1+distRatioCmp]);

    ylim([vert_top vert_bottom]); % 1 hPa is OK
    yticks(flip(vert_bottom:-100:vert_top));
    xticks(linspace(-90, 90, 7));
    
    % Set font properties...
    set(gca, 'FontName', 'Arial');
    set(gca, 'FontSize', 16);
    set(findall(gcf,'type','text'), 'FontName', 'Arial');
    set(findall(gcf,'type','text'), 'FontSize', 16);
    
    
end

sgtitle(sprintf(title_nice, var_nice_name, third_idx, fourth_idx));

% Set font properties...
set(gca, 'FontName', 'Arial');
set(gca, 'FontSize', 16);
set(findall(gcf,'type','text'), 'FontName', 'Arial');
set(findall(gcf,'type','text'), 'FontSize', 16);

% Set printing properties
%set(gcf, 'PaperUnits', 'inches', 'PaperPosition', [0 0 12 10]);
%set(gcf, 'Renderer', 'painters', 'Position', [90 90 1800 500])

set(gcf, 'Renderer', 'painters', 'Position', [600 800 2100 580]);
