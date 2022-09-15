%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%            _____            ______        %
% _______ _____  /_______________  /_______ %
% __  __ `__ \  __/  __ \  __ \_  /__  ___/ %
% _  / / / / / /_ / /_/ / /_/ /  / _(__  )  %
% /_/ /_/ /_/\__/ \____/\____//_/  /____/   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       "mtools" Research Toolkit           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Metrics_GlobalMass.m
%
% Compute global mass for an arbitrary set of data,
% given the coordinates and list of species.
% Also plots O3 column against OMI/MLS climatology
% if the O3 species is available.
%
% (c) 2019-2021 Haipeng Lin <jimmie.lin@gmail.com>
%
% Version: 2021.10.02
% Started: 2021.09.26
% See changelog at end of script
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% addpath("/usr/local/MATLAB/R2021a/toolbox/m_map");

%% data source configuration
% specify data coordinates. provide the "load" handle to the rest of the
% routines
coords_CESM = load('mdl_coords_CESM_f19_f19_mg17.mat');
coords_left = coords_CESM;
coords_right= coords_CESM;

fname_left  = "mdl_cesm2.1gc13.4.1_2016_monthly.mat";
fname_right = "mdl_cesm2.1gc13.4.1_2016_monthly.mat";

file_left   = load(fname_left);
file_right  = load(fname_right);

left_name      = "CESM2.1-GC13.4.1";
right_name     = "CESM2-GC";

%% data processing configuration (separated into section for ease of use)
% specify data indexing properties
%
% fourth_idx: usually date slice.
% if this is monthly data (indexed 1:12), set to -1 for a yearly average.
fourth_idx = 7; % date index

% specify minimum / maximum pressure properties
%         set ---                 minimum prs.        maximum prs. [Pa]
% for 100 hPa and above:            -999                 10 000
% for below 100 hPa    :           10 000               999 999
% for whole atm        :           
minPrsAtEdge = -999;
maxPrsAtEdge = 999999;

% species to compare.
spec_list = ["NOx", "NOy", "OH", "CH2O", "CO", "BrO", "ClO", "Cl", ...
              "HCl", "HBr", "HOCl", "HOBr", "CH3Cl", "BrCl", "BrNO3", "ClNO3", ...
              "HO2", "H2O2", "CH4", "PM25", "O3", ...
              "SO2", ...
              "NO", "NO2", "PAN", "HNO3", "NO3", "N2O", "N2O5", "HNO4", ...
              "MOH", "EOH", "ALD2", "C3H8", "DMS", "ACET", "MEK", "MVK", "TOLU", "MACR", "ALK4", "RCHO", "ISOP", ...
              "BCPI", "BCPO", "SO4", "OCPI", "OCPO"];
%spec_list = ["O3"];

%%%%%%%%%%%%%%%%%%%%%% NO USER CONFIGURABLE CODE BELOW %%%%%%%%%%%%%%%%%%%%
% read the coordinate data into separate vars, as necessary
lons_left   = coords_left.lons;
lats_left   = coords_left.lats;

lons_right  = coords_right.lons;
lats_right  = coords_right.lats;

% error check
if ~isfield(coords_left, 'hyai')
    error("ERROR: There is no hyai variable in the coordinates (left) struct. For this script to work, hyai coordinates MUST be provided.");
end

if ~isfield(coords_right, 'hyai')
    error("ERROR: There is no hyai variable in the coordinates (right) struct. For this script to work, hyai coordinates MUST be provided.");
end

if ~isfield(file_left, 'PSFC')
    error("ERROR: There is no PSFC variable in the data (left) struct. For this script to work, sfc pressure MUST be provided.");
end

if ~isfield(file_right, 'PSFC')
    error("ERROR: There is no PSFC variable in the data (right) struct. For this script to work, sfc pressure MUST be provided.");
end

% sort the species list in order
spec_list = sort(spec_list);

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

area_m2_left = area_m2_calc(lons_left, lats_left);
area_m2_right = area_m2_calc(lons_right, lats_right);

% Test Area_M2 if necessary
% for CESM processing - if lons are at 0.0 to 360.0, need to transpose data
% into -180.0, 180.0. check
% if max(lons_left) > 180.0
%    area_m2_left = circshift(area_m2_left, ...
%                  size(lons_left, 1)/2, ...
%                  1); % data needs to be in (lat, lon) 
% 
%    lons_left = lons_left - 180.0;
% 
%    fprintf("cesm-circshift: shifted data_left, lons_left\n");
% end
% 
% axg = subplot_tight(1,2,1);
% m_proj('miller', 'long',[-180.0 180.0], 'lat', [-87.5 87.5]);
% 
% m_pcolor(lons_left, lats_left, area_m2_left);
% colorbar;
% 
% m_plotbndry('/usr/local/MATLAB/R2021a/toolbox/boundary/world', ...
%     'color', 'k', 'linewidth', 1.5) % range lon is -180.0 180.0 so correct
% m_grid;
% 
% t = title("Area M2 test (left) [m^2]");
% 
% axg = subplot_tight(1,2,2);
% m_proj('miller', ...
%        'long',[-180.0 180.0], ...
%        'lat', [-87.5 87.5]);
% 
% m_pcolor(lons_right, lats_right, area_m2_right);
% colorbar;
% 
% m_plotbndry('/usr/local/MATLAB/R2021a/toolbox/boundary/world', ...
%     'color', 'k', 'linewidth', 1.5) % range lon is -180.0 180.0 so correct
% m_grid;
% 
% t = title("Area M2 test (right) [m^2]");

% convert MMR to VMR
% VMR = 28.9644 / molar_mass * MMR
% MMR = molar_mass / 28.9644 * VMR
% a hack for a molar mass table ...
molarTable = containers.Map;
molarTable('ACET') = 58.09;
molarTable('ALD2') = 44.06;
molarTable('ALK4') = 58.12;
molarTable('BrO') = 95.904;
molarTable('BrCl') = 115.45;
molarTable('BrNO3') = 141.91;
molarTable('C3H8') = 44.11;
molarTable('CH2O') = 30.031;
molarTable('CH3Cl') = 50.49;
molarTable('CH4') = 16.04;
molarTable('Cl') = 35.45;
molarTable('ClNO3') = 97.45;
molarTable('ClO') = 51.4521;
molarTable('CO') = 28.01;
molarTable('DMS') = 62.13;
molarTable('EOH') = 46.07;
molarTable('H2O2') = 34.02;
molarTable('HBr') = 80.91;
molarTable('HCl') = 36.458;
molarTable('HNO3') = 63.01;
molarTable('HNO4') = 79.01;
molarTable('HO2') = 33.01;
molarTable('HOBr') = 96.91;
molarTable('HOCl') = 52.49;
molarTable('ISOP') = 68.12;
molarTable('MACR') = 70.10;
molarTable('MEK') = 72.11;
molarTable('MOH') = 32.05;
molarTable('MVK') = 70.09;
molarTable('N2O') = 44.013;
molarTable('N2O5') = 108.02;
molarTable('NO') = 30.01;
molarTable('NO2') = 46.0055;
molarTable('NO3') = 62.01;
molarTable('O3') = 48.00;
molarTable('OH') = 17.008;
molarTable('PAN') = 121.05;
molarTable('RCHO') = 58.09;
molarTable('SO2') = 64.066;
molarTable('SO4') = 96.06;
molarTable('TOLU') = 92.15;
molarTable('XYLE') = 106.18;

% dummies:
molarTable('BCPI') = 12.01;
molarTable('BCPO') = 12.01;
molarTable('OCPI') = 12.01;
molarTable('OCPO') = 12.01;
molarTable('NOx') = 1.00;
molarTable('NOy') = 1.00;
molarTable('PM25') = 1.00;

% species specific code below. do loop
% sum species by multiplying [kg] in grid box by MMR
LM_left  = size(coords_left.hyam, 1);
LM_right = size(coords_right.hyam, 1);

% get surface pressures...
if fourth_idx > 0
    psfc_left = file_left.PSFC(:,:,fourth_idx);
    psfc_right = file_right.PSFC(:,:,fourth_idx);
    
    % compute PEDGE
    pedge_left = pedge_calc(coords_left.hyai, coords_left.hybi, psfc_left);
    pedge_right = pedge_calc(coords_right.hyai, coords_right.hybi, psfc_right);

    % compute DELP_DRY, which is just dry pedge deltas
    delp_dry_left = delp_dry_calc(pedge_left);
    delp_dry_right = delp_dry_calc(pedge_right);

    % compute grid box mass for each ... using area_m2, delp_dry
    AD_left = ad_calc(delp_dry_left, area_m2_left);
    AD_right = ad_calc(delp_dry_right, area_m2_right);
elseif fourth_idx == -1
    % if computing yearly average, all data needs to be computed for all 12
    % months.
    IM_left = size(file_left.PSFC, 1);
    JM_left = size(file_right.PSFC, 2);
    
    IM_right = size(file_left.PSFC, 1);
    JM_right = size(file_right.PSFC, 2);
    
    pedge_left = zeros(IM_left, JM_left, LM_left+1, 12);
    pedge_right = zeros(IM_right, JM_right, LM_right+1, 12);
    
    delp_dry_left = zeros(IM_left, JM_left, LM_left, 12);
    delp_dry_right = zeros(IM_right, JM_right, LM_right, 12);
    
    AD_left = zeros(IM_left, JM_left, LM_left, 12);
    AD_right = zeros(IM_right, JM_right, LM_right, 12);
    
    % for each, compute PEDGE
    for mo = 1:12
       pedge_left(:,:,:,mo) = pedge_calc(coords_left.hyai, coords_left.hybi, file_left.PSFC(:,:,mo));
       pedge_right(:,:,:,mo) = pedge_calc(coords_right.hyai, coords_right.hybi, file_right.PSFC(:,:,mo));
       
       delp_dry_left(:,:,:,mo) = delp_dry_calc(pedge_left(:,:,:,mo));
       delp_dry_right(:,:,:,mo) = delp_dry_calc(pedge_right(:,:,:,mo));
       
       AD_left(:,:,:,mo) = ad_calc(delp_dry_left(:,:,:,mo), area_m2_left);
       AD_right(:,:,:,mo) = ad_calc(delp_dry_right(:,:,:,mo), area_m2_right);
    end
end
% if max(lons_left) > 180.0
%    psfc_left = circshift(psfc_left, ...
%                  size(lons_left, 1)/2, ...
%                  1); % data needs to be in (lat, lon) 
% 
%    lons_left = lons_left - 180.0;
% 
%    fprintf("cesm-circshift: shifted data_left, lons_left\n");
% end
% m_proj('miller', 'long',[-180.0 180.0], 'lat', [-87.5 87.5]);
% m_pcolor(lons_left, lats_left, psfc_left);
% colorbar;
% m_plotbndry('/usr/local/MATLAB/R2021a/toolbox/boundary/world', ...
%     'color', 'k', 'linewidth', 1.5) % range lon is -180.0 180.0 so correct
% m_grid;

fprintf("Global Mass: Time slice %d\n", fourth_idx);
fprintf("Total Mass [Gg]        %-10s               %-10s\n", left_name, right_name);
fprintf("Species                Left_Model               Right_Model             %s diff (L-R)\n", '%')
fprintf("---------------------------------------------------------------------------------------\n")
for idx = 1:size(spec_list, 2)
    % zero
    Msum_l   = 0.0;
    Msum_r   = 0.0;

    % spec_current_left, spec_current_right
    spec_current_left  = file_left.(spec_list(idx));
    spec_current_right = file_right.(spec_list(idx));
    
    % spec_list(idx)
    MW = molarTable(spec_list(idx));
    VMRtoMMR = MW / 28.9644;
    
    for L = 1:LM_left
    %for L = 1:LM_left
        % sum by levels ..                 vv note this is element-wise mult.
        if fourth_idx > 0
            by_layer = sum(sum(AD_left(:,:,L) .* spec_current_left(:,:,L,fourth_idx) * VMRtoMMR));
        elseif fourth_idx == -1 % assuming fourth_idx spans from 1:12 ... yearly sum
            by_layer = 0; % accum.
            for mo = 1:12
                by_layer = by_layer + sum(sum(AD_left(:,:,L,mo) .* spec_current_left(:,:,L,mo) / 12 * VMRtoMMR));
            end
        end
        Msum_l = Msum_l + by_layer;
        %fprintf("=> [L] spc %s, L = %d, AD = %4.4e, v/v = %8.6e, m = %8.6f Gg\n", spec_list(idx), L, AD_left(25,25,L), spec_current_left(25,25,L,fourth_idx), by_layer/1e6);
    end
   
    for L = 1:LM_right
    %for L = 1:LM_right
        % sum by levels ..                  vv note this is element-wise mult.
        if fourth_idx > 0
            by_layer = sum(sum(AD_right(:,:,L) .* spec_current_right(:,:,L,fourth_idx) * VMRtoMMR));
        elseif fourth_idx == -1 % assuming fourth_idx spans from 1:12 ... yearly sum
            by_layer = 0; % accum.
            for mo = 1:12
                by_layer = by_layer + sum(sum(AD_right(:,:,L,mo) .* spec_current_right(:,:,L,mo) / 12 * VMRtoMMR));
            end
        end
        Msum_r = Msum_r + by_layer;
        %fprintf("=> [R] spc %s, L = %d, v/v = %8.6e, m = %8.6f Gg\n", spec_list(idx), L, vv, by_layer/1e6);
    end

    % fprintf("Column sums L = %.6f Gg, R = %.6f Gg\n", Msum_l/1e6, Msum_r/1e6)
    spc_name_nice = spec_list(idx);
    if MW == 1.0 || MW == 12.01
        spc_name_nice = append(spc_name_nice, " (dummy)");
    end
    fprintf("%-12s           %-8.6f \t\t%-8.6f \t\t%+02.3f\n", spc_name_nice, Msum_l/1e6, Msum_r/1e6, (Msum_l-Msum_r)/Msum_r*100)
end

% calculate ozone column value and plot (against climatology?)
% this is an added bonus IF the O3 field is present
% remember to put the processed omi_mls column ozone mat files in the same
% directory
%
% omi_mls_trop_column_ozone_5x5.mat
if isfield(file_left, 'O3') && isfield(file_right, 'O3') && false
    fprintf("=> Also showing ozone column plot...\n");
    
    omi_trop = load("omi_mls_trop_column_ozone_5x5.mat");
    omi_strat = load("omi_mls_strat_column_ozone_5x5.mat");
    % match fourth_idx ...
    if fourth_idx > 0
        omi_trop_data = omi_trop.out_mat(:,:,fourth_idx);
        omi_strat_data = omi_strat.out_mat(:,:,fourth_idx);
    
        omi_trop_data(omi_trop_data < 0) = NaN;
        omi_strat_data(omi_strat_data < 0) = NaN;
        
        omi_total_data = omi_trop_data + omi_strat_data;
    elseif fourth_idx == -1
        % yearly avg. 
        omi_trop_data = sum(omi_trop.out_mat(:,:,1:12), 3) / 12;
        omi_strat_data = sum(omi_strat.out_mat(:,:,1:12), 3) / 12;
    
        omi_trop_data(omi_trop_data < 0) = NaN;
        omi_strat_data(omi_strat_data < 0) = NaN;
    
        omi_total_data = omi_trop_data + omi_strat_data;
    end
    omi_total_data(omi_total_data < 0) = NaN;
    
    spec_current_left  = file_left.O3;
    spec_current_right = file_right.O3;
    
    % 1 DU = 2.1414e-5 kg/m2
    % from kg/m2 to DU: 1/2.1414e-5 = 46698 DU per kg/m2
    % calculate kg/m2 totals...
    IM_left = size(spec_current_left, 1);
    JM_left = size(spec_current_left, 2);
    
    IM_right = size(spec_current_right, 1);
    JM_right = size(spec_current_right, 2);
    
    O3_col_l = zeros(IM_left,  JM_left);
    O3_col_r = zeros(IM_right, JM_right);
    
    for L = 1:LM_left
        if fourth_idx > 0
            O3_col_l = O3_col_l + AD_left(:,:,L) .* spec_current_left(:,:,L,fourth_idx) ./ area_m2_left(:,:) * 46698 * 48 / 28.9644;
        elseif fourth_idx == -1
            for mo = 1:12
                O3_col_l = O3_col_l + AD_left(:,:,L,mo) .* spec_current_left(:,:,L,mo) / 12 ./ area_m2_left(:,:) * 46698 * 48 / 28.9644;
            end
        end
    end
    
    for L = 1:LM_right
        if fourth_idx > 0
            O3_col_r = O3_col_r + AD_right(:,:,L) .* spec_current_right(:,:,L,fourth_idx) ./ area_m2_right(:,:) * 46698 * 48 / 28.9644;
        elseif fourth_idx == -1
            for mo = 1:12
                O3_col_r = O3_col_r + AD_right(:,:,L,mo) .* spec_current_right(:,:,L,mo) / 12 ./ area_m2_right(:,:) * 46698 * 48 / 28.9644;
            end
        end
    end
    
    
    % sneakily auto fix cesm-circshift for O3_col_l and O3_col_r
    % handling of coordinates ...
    % cesm-circshift
    if max(lons_left) > 180.0
       O3_col_l = circshift(O3_col_l, ...
                     size(lons_left, 1)/2, ...
                     1); % data needs to be in (lat, lon) 

       lons_left_sf = lons_left - 180.0;

       fprintf("cesm-circshift: shifted data_left, lons_left\n");
    else
       lons_left_sf = lons_left;
    end

    if max(lons_right) > 180.0
       O3_col_r = circshift(O3_col_r, ...
                     size(lons_right, 1)/2, ...
                     1); % data needs to be in (lat, lon) 

       lons_right_sf = lons_right - 180.0;

       fprintf("cesm-circshift: shifted data_right, lons_right\n");
    else
       lons_right_sf = lons_right;
    end
    
    % plot column
    % LEFT
    axg = subplot_tight(2,3,1);

    hold on;
    m_proj('miller', ...
           'long',[-180.0 180.0], ...
           'lat', [-87.5 87.5]);

    m_pcolor(lons_left_sf, lats_left, O3_col_l);
    colorbar('southoutside');
    caxis([150, 450]);

    m_plotbndry('/usr/local/MATLAB/R2021a/toolbox/boundary/world', ...
        'color', 'k', 'linewidth', 1.5)
    m_grid;
    colormap(axg, flip(othercolor('Spectral6', 12)));

    t = title(sprintf(left_name));

    % RIGHT
    axg = subplot_tight(2,3,2);

    hold on;
    m_proj('miller', ...
           'long',[-180.0 180.0], ...
           'lat', [-87.5 87.5]);

    m_pcolor(lons_right_sf, lats_right, O3_col_r);
    colorbar('southoutside');
    caxis([150, 450]);

    m_plotbndry('/usr/local/MATLAB/R2021a/toolbox/boundary/world', ...
        'color', 'k', 'linewidth', 1.5)
    m_grid;
    colormap(axg, flip(othercolor('Spectral6', 12)));

    t = title(sprintf(right_name));
    
    
    % OMI/LMS 2004-2010 climatology
    axg = subplot_tight(2,3,3);

    hold on;
    m_proj('miller', ...
           'long',[-180.0 180.0], ...
           'lat', [-87.5 87.5]);

    m_pcolor(omi_trop.lons, omi_trop.lats, omi_total_data);
    colorbar('southoutside');
    caxis([150, 450]);

    m_plotbndry('/usr/local/MATLAB/R2021a/toolbox/boundary/world', ...
        'color', 'k', 'linewidth', 1.5)
    m_grid;
    colormap(axg, flip(othercolor('Spectral6', 12)));

    t = title("OMI/MLS 2004-2010 Climo (Total)");
    
    % do deltas with OMI/MLS climo. first do regridding of data to the
    % coarse OMI/MLS grid
    data_l2o = griddata(lons_left_sf, lats_left, O3_col_l, ...
                        omi_trop.lons, omi_trop.lats);
    data_r2o = griddata(lons_right_sf, lats_right, O3_col_r, ...
                        omi_trop.lons, omi_trop.lats);
    delta_l2o = data_l2o - omi_total_data;
    delta_r2o = data_r2o - omi_total_data;
    
    axg = subplot_tight(2,3,4);

    hold on;
    m_proj('miller', ...
           'long',[-180.0 180.0], ...
           'lat', [-87.5 87.5]);

    m_pcolor(omi_trop.lons, omi_trop.lats, delta_l2o);
    colorbar('southoutside');
    caxis([-45, 45]);

    m_plotbndry('/usr/local/MATLAB/R2021a/toolbox/boundary/world', ...
        'color', 'k', 'linewidth', 1.5)
    m_grid;
    colormap(axg, flip(othercolor('RdBu11', 9)));

    t = title(sprintf("%s - OMI/MLS Climo", left_name));
    
    
    axg = subplot_tight(2,3,5);

    hold on;
    m_proj('miller', ...
           'long',[-180.0 180.0], ...
           'lat', [-87.5 87.5]);

    m_pcolor(omi_trop.lons, omi_trop.lats, delta_r2o);
    colorbar('southoutside');
    caxis([-45, 45]);

    m_plotbndry('/usr/local/MATLAB/R2021a/toolbox/boundary/world', ...
        'color', 'k', 'linewidth', 1.5)
    m_grid;
    colormap(axg, flip(othercolor('RdBu11', 9)));

    t = title(sprintf("%s - OMI/MLS Climo", right_name));
    
    sgtitle(sprintf("O3 column [DU] - date slice %d", fourth_idx));

    set(gcf, 'Renderer', 'painters', 'Position', [90 0 1250 970]);
end

%%%%%%%%%%%%%%%%%%%%%% GENERIC UTILITIES CODE BELOW %%%%%%%%%%%%%%%%%%%%

% Routine utilities for computing air quantities
% ported from the WRF-GC coupler, GEOS-Chem gc_grid_mod.F90,
% GEOS-Chem calc_met_mod.F90, among others
% (hplin, 9/26/21)

% AD: calculates air mass for each level
% dry air mass [kg]
function AD = ad_calc(DELP_DRY, AREA_M2)
    % from physconstants.F90 (G-C)
    G0 = 9.80665;
    % G0_100 = 100.0 / G0;
    
    IM = size(DELP_DRY, 1);
    JM = size(DELP_DRY, 2);
    LM = size(DELP_DRY, 3);
    
    AD = zeros(IM, JM, LM);
    for L = 1:LM
        %            note element-wise mult vv
        AD(:,:,L) = DELP_DRY(:,:,L) .* AREA_M2(:,:) / G0; % [Pa] calc.
        % if using [hPa], need to use G0_100 like in G-C
    end
end

% pedge: calculates EDGE pressure levels given hyai, hybi, PSFC
% for the entire grid.
function pedge = pedge_calc(hyai, hybi, PSFC)
    % sizes
    IM = size(PSFC, 1);
    JM = size(PSFC, 2);
    LM = size(hyai, 1) - 1;
    %DT = size(PSFC, 3);
    
    % assumes PSFC is in [Pa] (but outputs just follow PSFC units anyway)
    % if you want to use this for delp_dry calc and AD (kg) calc though
    % must give Pa
    % pedge = zeros(IM, JM, LM+1, DT);
    pedge = zeros(IM, JM, LM+1);
    
    fprintf("\n=> Computing PEDGE ... IM=%d,JM=%d,LM=%d\n",IM,JM,LM)
    
    %for T = 1:DT
    for L = 1:LM+1
        % pedge(:,:,L,T) = hyai(L) + hybi(L) .* PSFC(:,:,T);
        
        % GEOS-Chem format
        %pedge(:,:,L) = hyai(L) + (hybi(L) * PSFC(:,:));
        
        % CAM-chem format https://www.ncl.ucar.edu/Document/Functions/Built-in/pres_hybrid_ccm.shtml
        pedge(:,:,L) = hyai(L) * 100000.0 + hybi(L) * PSFC(:,:);
    end
    %end
end

% delp_dry: give PEDGE_DRY, returns deltas for each grid box.
function delp_dry = delp_dry_calc(PEDGE)
    IM = size(PEDGE, 1);
    JM = size(PEDGE, 2);
    LM = size(PEDGE, 3) - 1;
    
    delp_dry = zeros(IM, JM, LM);
    for L = 1:LM
        delp_dry(:,:,L) = PEDGE(:,:,L) - PEDGE(:,:,L+1);
    end
end

% area_m2: calculates M2 area of grid box given CENTER lons, lats
% which have been expanded to 2-D sizes lons=lats=(IM,JM)
%
% assumes CENTER because that is what we work with out of COARDS/CF
% conventions. cf. https://cfconventions.org/cf-conventions/cf-conventions.html
% "Bounds for 2-D coordinate variables with 4-sided cells"
%
% currently only works on cartesian grid.
%
% assumes lats are from south to north
function area_m2 = area_m2_calc(lons, lats)
    % Re: Radius of Earth [m]
    Re = 6.3710072e6;
    
    % assuming these are center quantities.
    IM = size(lons, 1);
    JM = size(lats, 2);
    
    % create area_m2:
    % fixme: may not be exactly accurate for curvilinear 
    % due to the calculation of dx is
    % based on adjacent center cells.
    area_m2 = zeros(IM, JM);
    
    % loop over the grid
    for J = 1:JM
        for I = 1:IM
            % compute the edge radian of the y-coordinate to north and
            % south. given the centers, we can do 1/2
            %if lons(I,J) == 90.0
            %    YEdgeRad_N = pi/2;
            %elseif lons(I,J) == -90.0
            %    YEdgeRad_N = -pi/2;
            %else
           if J == 1
               YEdge_N = (lats(I,J+1) + lats(I,J))/2;
               YEdge_S = -90.0;
           elseif J == JM
               YEdge_N =  90.0;
               YEdge_S = (lats(I,J-1) + lats(I,J))/2;
           else
               YEdge_N = (lats(I,J+1) + lats(I,J))/2;
               YEdge_S = (lats(I,J-1) + lats(I,J))/2;
           end
           
           YEdgeRad_N = YEdge_N * pi / 180;
           YEdgeRad_S = YEdge_S * pi / 180;
           
           YEdgeSin_N = sin(YEdgeRad_N);
           YEdgeSin_S = sin(YEdgeRad_S);
           
           % get adjacent grid box dx (assuming rectilinearity here...)
           if I ~= IM
               DX = lons(I+1,J) - lons(I,J);
           else
               DX = lons(I,J)   - lons(I-1,J);
           end
           
           area_m2(I,J) = (DX * pi / 180) * (Re^2) * (YEdgeSin_N - YEdgeSin_S);
        end
    end
    

end
