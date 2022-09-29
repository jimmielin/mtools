%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%            _____            ______        %
% _______ _____  /_______________  /_______ %
% __  __ `__ \  __/  __ \  __ \_  /__  ___/ %
% _  / / / / / /_ / /_/ / /_/ /  / _(__  )  %
% /_/ /_/ /_/\__/ \____/\____//_/  /____/   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       "mtools" Research Toolkit           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot_MDLvOBS.m
%
% Plot model outputs against observations (ground)
% (c) 2019-2022 Haipeng Lin <work@jimmielin.me>
%
% Version: 2022.09.28
% See changelog at end of script
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Configuration options
% Metadata - will be used to get species data (magic) as well
% meta_title: %d inserts date slice, %s inserts species name
var_name     = "O3";
meta_title   = "2016 monthly average - date slice %d - species %s [China MEP]";

% Input data files
mdl_filename = 'mdl_cesm2.1gc13.4.1_2016_monthly.mat';
mdl_CESM2GC = load(mdl_filename, var_name);
coords_CESM_f19_f19_mg17 = load("mdl_coords_CESM_f19_f19_mg17.mat");

obs_chn     = load("obs_chinamep_ground/2016_hourly_all.mat");
obs_kor     = load("obs_airkorea_ground/2016_hourly_all.mat");

% figure output location (no trailing underscore needed)
figure_out_prefix = "out_figures/cesm2.1gc13.4.1_2016_chinamep";

% specify data indexing properties
%
% third_idx:  ignored in this script since only ground obs.
% fourth_idx: usually date slice
%
% WARNING: IF DATA IS GIVEN IN 3-DIMENSIONAL ARRAY, USUALLY FOR 2-D DATA
% IT WILL BE INTERPRETED AS (lon, lat, fourth_idx)
% BECAUSE IT USUALLY MEANS IT IS 2-D DATA, AND FOURTH_IDX IS DATE SLICE BY
% CONVENTION.
third_idx =  1; % level index (ignored here)
fourth_idx = 7; % date index

% specify observation data indexing properties, if computing averages
% usually, obs_month is fourth_idx.

% obs_mode: averaged
obs_mode    = "averaged";

obs_current = obs_chn;
obs_year    = 2016;
obs_month   = fourth_idx;

% Read data using magic
mdl_coords   = coords_CESM_f19_f19_mg17;
eval(sprintf("mdl_data     = mdl_CESM2GC.%s;", var_name));
eval(sprintf("obs_data     = obs_current.obs_%s;", lower(var_name)));

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
mdl_handling = "cesm-circshift";

% specify data scaling properties (generally does not need change)
% in CESM, VMR is in 1 so convert to ppm in 1e6
convert_factor = 1e6;

% Mode:
% - closest_model_for_obs - for each observation point find closest model
% grid box.
% - average_obs_over_model - average observational points over model
compare_mode = "average_obs_over_model";

%% No user configurable options below
% month end
obs_month_end = 31;
if obs_month == 2
    obs_month_end = 28;
    if (mod(obs_year, 4) == 0 && mod(obs_year, 100) ~= 0) || mod(obs_year, 400) == 0
        obs_month_end = 29;
    end
elseif obs_month == 4 || obs_month == 6 || obs_month == 9 || obs_month == 11
    obs_month_end = 30;
end

% handling of coordinates ...
mdl_coords_lon = mdl_coords.lons;
mdl_coords_lat = mdl_coords.lats;
mdl_coords_lev = mdl_coords.levs;

% cesm-circshift
if strcmp(mdl_handling, 'cesm-circshift')
    % for CESM processing - if lons are at 0.0 to 360.0, need to transpose data
    % into -180.0, 180.0. check
    if max(mdl_coords_lon) > 180.0
       mdl_data = circshift(mdl_data, ...
                     size(mdl_coords_lon, 1)/2, ...
                     1); % data needs to be in (lat, lon) 

       mdl_coords_lon = mdl_coords_lon - 180.0;

       fprintf("cesm-circshift: shifted mdl_data, mdl_coords_lon\n");
    end
end

% auto-fix 2-D coordinates
if size(mdl_coords_lon, 2) == 1
    lons_new = zeros(size(mdl_coords_lon, 1), size(mdl_coords_lat, 1));
    lats_new = zeros(size(mdl_coords_lon, 1), size(mdl_coords_lat, 1));

    for j = 1:size(mdl_coords_lat, 1)
        lons_new(:,j) = mdl_coords_lon(:,1);
    end

    for i = 1:size(mdl_coords_lon, 1)
        lats_new(i,:) = mdl_coords_lat(:,1);
    end

    mdl_coords_lon = lons_new;
    mdl_coords_lat = lats_new;
end

% Slice the data subset
if length(size(mdl_data)) == 4
    if third_idx > 0
        mdl_data  = mdl_data (:,:,third_idx,fourth_idx) * convert_factor;
    else
        % fancy processing - if third_idx < 1, collapse the third_idx
        mdl_data  = sum(mdl_data (:,:,:,fourth_idx), 3) * convert_factor;
    end
elseif length(size(mdl_data)) == 3
    fprintf("Note: 2-D data detected; data is interpreted as lon,lat,fourth_idx. buyer beware.\n")
    
    mdl_data  = mdl_data (:,:,fourth_idx) * convert_factor;
else
    error("size of mdl_data is neither 3-D or 4-D. what is going on?");
end

% create averaged observations if needed, for monthly
if obs_mode == "averaged"
    % build monthly average. to make more generic later
    % find beginning and end time indices
    obs_time_slice_start = obs_current.date2ptr(sprintf("%04d%02d0101", obs_year, obs_month));
    obs_time_slice_end   = obs_current.date2ptr(sprintf("%04d%02d%02d24", obs_year, obs_month, obs_month_end));
    obs_data_slice = obs_data(:, obs_time_slice_start:obs_time_slice_end);
    obs_data = nanmean(obs_data_slice, 2);
end

% Now data is subsliced into mdl_data(i, j)
serial_set_obs = [];
serial_set_mdl = [];

if compare_mode == "closest_model_for_obs"
    % Compare against observations. Loop for each observational data point...
    
    % expects that data in obs_data is now 1d
    for obs_siteidx = 1:1:obs_current.sitenum
        site  = obs_current.sites(obs_siteidx, 1);
        slong = obs_current.sites(obs_siteidx, 2);
        slat  = obs_current.sites(obs_siteidx, 3);
        
        % prelim filter out stations with NaN data
        if isnan(obs_data(obs_siteidx, 1))
            continue
        end
        %fprintf("Site IntID %d, SiteID %6d, Long x Lat %f x %f\n", obs_siteidx, site, slong, slat);
    
        % This assumes rectilinearity in xlat, xlong, also not very efficient as
        % it goes to the end of the loop
        ddmax = 360;
        di = 1;
        dj = 1;
        for i = 1:size(mdl_coords_lon, 1)
            if abs(mdl_coords_lon(i,1)-slong) < ddmax
                ddmax = abs(mdl_coords_lon(i,1)-slong);
                di = i;
            end
        end
    
        ddmax = 180;
        for j = 1:size(mdl_coords_lat, 2)
            if abs(mdl_coords_lat(1,j)-slat) < ddmax
                ddmax = abs(mdl_coords_lat(1,j)-slat);
                dj = j;
            end
        end
    
        % disallow ddmax by threshold because if it is too far it has no
        % meaning (0.01 deg ~ 1 km)
        if ddmax > 0.1
            continue
        end
    
        %fprintf("=> Closest grid box i,j,mdl_long,mdl_lat = %d, %d, %f, %f\n", di, dj, mdl_coords_lon(di,dj), mdl_coords_lat(di,dj));
    
        % Now that the site is pinpointed, compare against the corresponding
        % mdl_data(i,j), which has been converted to ug/m3 (pm) or ppm (others)
        serial_set_obs = [serial_set_obs; obs_data(obs_siteidx, 1)];
        serial_set_mdl = [serial_set_mdl; mdl_data(di, dj)];
      
    end
elseif compare_mode == "average_obs_over_model"
    % loop for each grid box...
    % measure deltas. assume some rectilinearity here given these are
    % lat/lon centers
    dlon = abs(mdl_coords_lon(2,1) - mdl_coords_lon(1,1));
    dlat = abs(mdl_coords_lat(1,2) - mdl_coords_lat(1,1));
    for i = 1:size(mdl_coords_lon, 1)
        for j = 1:size(mdl_coords_lat, 2)
            % find obs to average... this could be more efficient
            averaged_obs = [];
            for obs_siteidx = 1:1:obs_current.sitenum
                site  = obs_current.sites(obs_siteidx, 1);
                slong = obs_current.sites(obs_siteidx, 2);
                slat  = obs_current.sites(obs_siteidx, 3);

                % are we in grid box? check
                if slong >= mdl_coords_lon(i,j)-dlon && slong <=mdl_coords_lon(i,j)+dlon && ...
                        slat >=mdl_coords_lat(i,j)-dlat && slat<=mdl_coords_lat(i,j)+dlat
                    averaged_obs = [averaged_obs; obs_data(obs_siteidx, 1)];
                end
            end
            
            obs_count = size(averaged_obs, 1);
            if obs_count > 0
                averaged_obs = nanmean(averaged_obs);
                serial_set_obs = [serial_set_obs; averaged_obs];
                serial_set_mdl = [serial_set_mdl; mdl_data(i, j)];
                fprintf("=> found %d observations in %d, %d\n", obs_count, i, j);
            end
        end
    end
end

%% Make plot
scatter(serial_set_obs, serial_set_mdl, 16, 'filled');
hold on;
plot([0, 1], [0, 1], '--', 'LineWidth', 2, 'color', [0 0 0]+0.72);
xlim([0 max(serial_set_mdl)]);
ylim([0 max(serial_set_mdl)]);

coef_r = corrcoef(serial_set_obs, serial_set_mdl);
coef_r = coef_r(2);

% linreg. note, x has error. should fix this later
%X = [ones(length(serial_set_obs), 1) serial_set_obs];
%b = X\serial_set_mdl;
%yCalc2 = X*b;
%plot(serial_set_obs, yCalc2, '--');
%legend('data', '1:1', sprintf('slope = %.2f, intercept = %.2f\nr = %.2f', b(1), b(2), coef_r));

% generalized major axis regression: requires gmregress file
[b,bintr,bintjm] = gmregress(serial_set_obs, serial_set_mdl);
dummy_xs = linspace(0, 1, 100);
yCalc2 = dummy_xs * b(2) + b(1);
plot(dummy_xs, yCalc2, 'LineWidth', 2);
    % ±X
legend(sprintf('data (%d points)', length(serial_set_obs)), '1:1', ...
       sprintf('slope = %.2f, intercept = %.2f\nr = %.2f', ...
               b(2), b(1), coef_r), 'Location', 'best');

title(sprintf(meta_title, fourth_idx, var_name), 'Interpreter', 'none');
subtitle(mdl_filename, 'Interpreter', 'none');
xlabel("Observations [ppm]");
ylabel("Model [ppm]");

set(gca, 'FontName', 'Arial');
set(gca, 'FontSize', 16);
set(findall(gcf,'type','text'), 'FontName', 'Arial');
set(findall(gcf,'type','text'), 'FontSize', 16);

set(gcf, 'Renderer', 'painters', 'Position', [300 300 800 800]);

fname = sprintf("%s_compare_%s_date%02d.png", figure_out_prefix, var_name, fourth_idx);
saveas(gcf, fname);
close(gcf);

fprintf("saved to %s\n", fname);

% Changelog:
% 2022.09.28 - Average observational points over model grid instead of
%              the other way around.
%            - Save plots to .png consistent with other mtools.
% 2022.09.13 - Initial version