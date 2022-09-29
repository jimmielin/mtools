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
% Also can read a set of reference mass values from
% Util_ReadReferenceGlobalMass.m if available,
% used to compare against GEOS-Chem benchmarking
% results in FullChem/BenchmarkResults/Tables/\
% GlobalMass_{Trop,TropStrat}_{Jan,Apr,Jul,Oct}2019.txt
%
% (c) 2019-2022 Haipeng Lin <jimmie.lin@gmail.com>
%
% Version: 2022.09.15
% Started: 2021.09.26
% See changelog at end of script
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath("/usr/local/MATLAB/R2021b/toolbox/m_map");

%% data source configuration
% specify data coordinates. provide the "load" handle to the rest of the
% routines
coords_f19 = load('mdl_coords_CESM_f19_f19_mg17.mat');
coords_f09 = load('mdl_coords_CESM_f09_f09_mg17.mat');
coords_left = coords_f19;
coords_right = coords_f19;
%coords_right= coords_f09;

fname_left  = "mdl_cesm2.1gc13.4.1_2016_monthly.mat";
%fname_right = "mdl_cesm2.2camchem_2016_monthly.mat";
fname_right = "mdl_cesm2.1gc13.1.2_2016_monthly.mat";

file_left   = load(fname_left);
file_right  = load(fname_right);

left_name      = "CESM2.1-GC13.4.1";
%right_name     = "CESM2.2-CAMChem";
right_name = "CESM2.1-GC13.1.2";

% now also have a reference model output
fname_reference = "mdl_ref_gc_benchmark/mdl_ref_gc_benchmark_13.4.0rc.4-tropstrat-201907.mat";
load(fname_reference, "mdl_ref_sum_map");

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

LM_left  = size(coords_left.hyam, 1);
LM_right = size(coords_right.hyam, 1);
LM_left_range   = 1:LM_left;
LM_right_range  = 1:LM_right;

% figure output location (no trailing underscore needed)
figure_out_prefix = sprintf("out_figures/%s_%s_2016", left_name, right_name);

% species to compare.
spec_list = ["NOx", "NOy", "OH", "CH2O", "CO", "BrO", "ClO", "Cl", ...
              "HCl", "HBr", "HOCl", "HOBr", "CH3Cl", "BrCl", ... %"BrNO3", "ClNO3", ...
              "HO2", "H2O2", "CH4", "PM25", "O3", ...
              "SO2", "SO4", ...
              "NO", "NO2", "PAN", "HNO3", "NO3", "N2O", "N2O5", ... %"HNO4", ...
              "MOH", "EOH", "ALD2", "C3H8", "DMS", "ACET", "MEK", "MVK", "TOLU", "MACR", "ALK4", "ISOP", ...
              "BCPI", "BCPO", "OCPI", "OCPO"];
% simplified old list:
spec_list = ["BrO", "CH2O", "CH3Cl", "CH4", "ClO", "CO", "HCl", "HO2", "HOCl", "ISOP", "N2O", "NO", "NO2", "O3", "OH", "PAN", "SO2"];

% omi_mls mode: all|trop|strat
omi_mls_mode = "all";

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
%
% note that CESM data is in VMR (mol/mol), not ppmv, so to ppm scale by 1e6
%
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

fprintf("Global Mass: Time slice %d \t\t LM_left: %d:%d (~%.2f to %.2f hPa) \t\t LM_right: %d:%d (~%.2f to %.2f hPa)\n", fourth_idx, ...
        LM_left_range(1), LM_left_range(end), mean(pedge_left(:,:,LM_left_range(1)),'all')/100, mean(pedge_left(:,:,LM_left_range(end))/100,'all'), ...
        LM_right_range(1), LM_right_range(end), mean(pedge_right(:,:,LM_right_range(1)),'all')/100, mean(pedge_right(:,:,LM_right_range(end))/100,'all'));
fprintf("Reference model: %s\n", fname_reference);
fprintf("Total Mass [Gg]        %-15s         %-15s\n", left_name, right_name);
fprintf("Species                Left_Model               Right_Model             %s diff (L-R)          Reference_Model\n", '%')
fprintf("--------------------------------------------------------------------------------------------------------------\n")
for idx = 1:size(spec_list, 2)
    % zero
    Msum_l   = 0.0;
    Msum_r   = 0.0;

    % check if magic is required. e.g., special handling is performed for
    % SO4 if not present, where it is replaced with so4_a1+so4_a2+so4_a3;
    %
    % fixme: should actually +H2SO4 as well
    is_fallback_left = 0;
    is_fallback_right = 0;
    if spec_list(idx) == 'SO4'
        if ~isfield(file_left, 'SO4') && isfield(file_left, 'so4_a1')
            spec_current_left  = file_left.so4_a1 + file_left.so4_a2 + file_left.so4_a3;
            is_fallback_left = 1;
        end

        if ~isfield(file_right, 'SO4') && isfield(file_right, 'so4_a1')
            spec_current_right  = file_right.so4_a1 + file_right.so4_a2 + file_right.so4_a3;
            is_fallback_right = 1;
        end
    end

    % spec_current_left, spec_current_right
    if ~is_fallback_left
        spec_current_left  = file_left.(spec_list(idx));
    end
    if ~is_fallback_right 
        spec_current_right = file_right.(spec_list(idx));
    end
    
    % spec_list(idx)
    MW = molarTable(spec_list(idx));
    VMRtoMMR = MW / 28.9644;
    
    for L = LM_left_range
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
   
    for L = LM_right_range
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
    
    pct_diff = (Msum_l-Msum_r)/Msum_r*100;
    extra = "   ";
    if abs(pct_diff) > 20
        extra = "!!!";
    elseif abs(pct_diff) > 10
        extra = "!! ";
    elseif abs(pct_diff) > 5
        extra = "! ";
    end

    reference = -1.0;
    if isKey(mdl_ref_sum_map, spec_list(idx))
        reference = mdl_ref_sum_map(spec_list(idx));
    end

    fprintf("%-12s           %-8.6f  \t\t%-8.6f \t\t%+8.3f %s\t\t%-8.6f\n", spc_name_nice, Msum_l/1e6, Msum_r/1e6, pct_diff, extra, reference)
end

%% calculate ozone column value and plot (against climatology?)
% this is an added bonus IF the O3 field is present
% remember to put the processed omi_mls column ozone mat files in the same
% directory
%
% omi_mls_trop_column_ozone_5x5.mat
if isfield(file_left, 'O3') && isfield(file_right, 'O3')
    fprintf("=> Also showing ozone column plot... mode = %s\n", omi_mls_mode);
    
    omi_trop = load("obs_omi_mls/omi_mls_climatology_trop_column_ozone_5x5.mat");
    omi_strat = load("obs_omi_mls/omi_mls_climatology_strat_column_ozone_5x5.mat");
    % match fourth_idx ...
    if fourth_idx > 0
        omi_trop_data = omi_trop.out_mat(:,:,fourth_idx);
        omi_strat_data = omi_strat.out_mat(:,:,fourth_idx);
    
        omi_trop_data(omi_trop_data < 0) = NaN;
        omi_strat_data(omi_strat_data < 0) = NaN;
        
        if omi_mls_mode == "trop"
            omi_total_data = omi_trop_data;
        elseif omi_mls_mode == "strat"
            omi_total_data = omi_strat_data;
        else
            omi_total_data = omi_trop_data + omi_strat_data;
        end
    elseif fourth_idx == -1
        % yearly avg. 
        omi_trop_data = sum(omi_trop.out_mat(:,:,1:12), 3) / 12;
        omi_strat_data = sum(omi_strat.out_mat(:,:,1:12), 3) / 12;
    
        omi_trop_data(omi_trop_data < 0) = NaN;
        omi_strat_data(omi_strat_data < 0) = NaN;
    
        if omi_mls_mode == "trop"
            omi_total_data = omi_trop_data;
        elseif omi_mls_mode == "strat"
            omi_total_data = omi_strat_data;
        else
            omi_total_data = omi_trop_data + omi_strat_data;
        end
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
    
    % hack hack
    if omi_mls_mode == "trop"
        LM_left_s = 1; LM_left_e = 35;
        LM_right_s = 1; LM_right_e = 35;

        omimls_o3_range = [15, 45];
        omimls_o3_delta_range = [-15, 15];
    elseif omi_mls_mode == "strat"
        LM_left_s = 36; LM_left_e = LM_left;
        LM_right_s = 36; LM_right_e = LM_right;

        omimls_o3_range = [150, 450];
        omimls_o3_delta_range = [-45, 45];
    else
        LM_left_s = 1; LM_left_e = LM_left;
        LM_right_s = 1; LM_right_e = LM_right;

        omimls_o3_range = [150, 450];
        omimls_o3_delta_range = [-45, 45];
    end

    for L = LM_left_s:LM_left_e
        if fourth_idx > 0
            O3_col_l = O3_col_l + AD_left(:,:,L) .* spec_current_left(:,:,L,fourth_idx) ./ area_m2_left(:,:) * 46698 * 48 / 28.9644;
        elseif fourth_idx == -1
            for mo = 1:12
                O3_col_l = O3_col_l + AD_left(:,:,L,mo) .* spec_current_left(:,:,L,mo) / 12 ./ area_m2_left(:,:) * 46698 * 48 / 28.9644;
            end
        end
    end
    
    for L = LM_right_s:LM_right_e
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
    caxis(omimls_o3_range);

    m_plotbndry('/usr/local/MATLAB/R2021b/toolbox/boundary/world', ...
        'color', 'k', 'linewidth', 1.5)
    m_grid;
    colormap(axg, flip(othercolor('Spectral6', 12)));

    t = title(sprintf(left_name));
    set(gca, 'FontName', 'Arial');
    set(gca, 'FontSize', 16);
    set(findall(gcf,'type','text'), 'FontName', 'Arial');
    set(findall(gcf,'type','text'), 'FontSize', 16);

    % RIGHT
    axg = subplot_tight(2,3,2);

    hold on;
    m_proj('miller', ...
           'long',[-180.0 180.0], ...
           'lat', [-87.5 87.5]);

    m_pcolor(lons_right_sf, lats_right, O3_col_r);
    colorbar('southoutside');
    caxis(omimls_o3_range);

    m_plotbndry('/usr/local/MATLAB/R2021b/toolbox/boundary/world', ...
        'color', 'k', 'linewidth', 1.5)
    m_grid;
    colormap(axg, flip(othercolor('Spectral6', 12)));

    t = title(sprintf(right_name));
    set(gca, 'FontName', 'Arial');
    set(gca, 'FontSize', 16);
    set(findall(gcf,'type','text'), 'FontName', 'Arial');
    set(findall(gcf,'type','text'), 'FontSize', 16);
    
    
    % OMI/LMS 2004-2010 climatology
    axg = subplot_tight(2,3,3);

    hold on;
    m_proj('miller', ...
           'long',[-180.0 180.0], ...
           'lat', [-87.5 87.5]);

    m_pcolor(omi_trop.lons, omi_trop.lats, omi_total_data);
    colorbar('southoutside');
    caxis(omimls_o3_range);

    m_plotbndry('/usr/local/MATLAB/R2021b/toolbox/boundary/world', ...
        'color', 'k', 'linewidth', 1.5)
    m_grid;
    colormap(axg, flip(othercolor('Spectral6', 12)));

    t = title("OMI/MLS 2004-2010 Climo (Total)");
    set(gca, 'FontName', 'Arial');
    set(gca, 'FontSize', 16);
    set(findall(gcf,'type','text'), 'FontName', 'Arial');
    set(findall(gcf,'type','text'), 'FontSize', 16);
    
    % do deltas with OMI/MLS climo. first do regridding of data to the
    % coarse OMI/MLS grid
    data_l2o = griddata(lons_left_sf, lats_left, O3_col_l, ...
                        omi_trop.lons, omi_trop.lats);
    data_r2o = griddata(lons_right_sf, lats_right, O3_col_r, ...
                        omi_trop.lons, omi_trop.lats);

    % also regrid left to right if necessary - check which is coarser
    if size(lons_left_sf, 1) > size(lons_right_sf, 1)
        % left is finer, regrid to right (coarse grid)
        data_r2c = O3_col_r;
        data_l2c = griddata(lons_left_sf, lats_left, O3_col_l, ...
                            lons_right_sf, lats_right);
        lons_c = lons_right_sf; lats_c = lats_right;
        fprintf("==> regridding left data to right (right is coarser)\n");
    elseif size(lons_left_sf, 1) < size(lons_right_sf, 1)
        % left is coarser, regrid to left
        data_l2c = O3_col_l;
        data_r2c = griddata(lons_right_sf, lats_right, O3_col_r, ...
                            lons_left_sf, lats_left);
        lons_c = lons_left_sf; lats_c = lats_left;
        fprintf("==> regridding R data to L (L is coarser)\n");
    else
        data_l2c = O3_col_l;
        data_r2c = O3_col_r;
        % both on same grid
        lons_c = lons_left_sf; lats_c = lats_left;
        fprintf("==> L and R are on same grid, no regridding.\n");
    end

    delta_l2o = data_l2o - omi_total_data;
    delta_r2o = data_r2o - omi_total_data;
    delta_lminusr = data_l2c - data_r2c;
    
    axg = subplot_tight(2,3,4);

    hold on;
    m_proj('miller', ...
           'long',[-180.0 180.0], ...
           'lat', [-87.5 87.5]);

    m_pcolor(omi_trop.lons, omi_trop.lats, delta_l2o);
    colorbar('southoutside');
    caxis(omimls_o3_delta_range);

    m_plotbndry('/usr/local/MATLAB/R2021b/toolbox/boundary/world', ...
        'color', 'k', 'linewidth', 1.5)
    m_grid;
    colormap(axg, flip(othercolor('RdBu11', 9)));

    t = title(sprintf("%s - OMI/MLS Climo", left_name));
    set(gca, 'FontName', 'Arial');
    set(gca, 'FontSize', 16);
    set(findall(gcf,'type','text'), 'FontName', 'Arial');
    set(findall(gcf,'type','text'), 'FontSize', 16);
    
    
    axg = subplot_tight(2,3,5);

    hold on;
    m_proj('miller', ...
           'long',[-180.0 180.0], ...
           'lat', [-87.5 87.5]);

    m_pcolor(omi_trop.lons, omi_trop.lats, delta_r2o);
    colorbar('southoutside');
    caxis(omimls_o3_delta_range);

    m_plotbndry('/usr/local/MATLAB/R2021b/toolbox/boundary/world', ...
        'color', 'k', 'linewidth', 1.5)
    m_grid;
    colormap(axg, flip(othercolor('RdBu11', 9)));

    t = title(sprintf("%s - OMI/MLS Climo", right_name));
    set(gca, 'FontName', 'Arial');
    set(gca, 'FontSize', 16);
    set(findall(gcf,'type','text'), 'FontName', 'Arial');
    set(findall(gcf,'type','text'), 'FontSize', 16);

    % Also do a diff of L-R while at it
    axg = subplot_tight(2,3,6);

    hold on;
    m_proj('miller', ...
           'long',[-180.0 180.0], ...
           'lat', [-87.5 87.5]);

    m_pcolor(lons_c, lats_c, delta_lminusr);
    colorbar('southoutside');
    caxis(omimls_o3_delta_range);

    m_plotbndry('/usr/local/MATLAB/R2021b/toolbox/boundary/world', ...
        'color', 'k', 'linewidth', 1.5)
    m_grid;
    colormap(axg, flip(othercolor('RdBu11', 9)));

    t = title(sprintf("%s - %s (Diff L-R)", left_name, right_name));
    
    sgtitle(sprintf("O3 column [DU] - date slice %d", fourth_idx));

    set(gcf, 'Renderer', 'painters', 'Position', [90 0 1600 1150]);
    set(gca, 'FontName', 'Arial');
    set(gca, 'FontSize', 16);
    set(findall(gcf,'type','text'), 'FontName', 'Arial');
    set(findall(gcf,'type','text'), 'FontSize', 16);
    
    fname = sprintf("%s_%02d_o3_omimls_column.png", figure_out_prefix, fourth_idx);
    saveas(gcf, fname);
    close(gcf);
end

%% calculate CO column value and plot (against MOPITT obs)
% this is an added bonus IF the O3 field is present
% remember to put the processed omi_mls column ozone mat files in the same
% directory
%
% omi_mls_trop_column_ozone_5x5.mat
if isfield(file_left, 'CO') && isfield(file_right, 'CO')
    fprintf("=> Also showing MOPITT CO column density plot...\n");

    mopitt_co_range = [0, 30e17];
    mopitt_co_delta_range = [-9e17, 9e17];
    
    mopitt_year = 2016;
    if fourth_idx > 0
        load(sprintf("obs_mopitt_co/%04d%02d_monthly_mean.mat", mopitt_year, fourth_idx));
    else
        load(sprintf("obs_mopitt_co/%04d_yearly_mean.mat", mopitt_year));
    end
    mopitt_co_column_day(mopitt_co_column_day < 0) = NaN;
    mopitt_co_column_night(mopitt_co_column_night < 0) = NaN;

    mopitt_co_column = transpose(mopitt_co_column_night); % huh? use day or night. fixme
    
    spec_current_left  = file_left.CO;
    spec_current_right = file_right.CO;
    
    % convert to molec/cm2
    IM_left = size(spec_current_left, 1);
    JM_left = size(spec_current_left, 2);
    
    IM_right = size(spec_current_right, 1);
    JM_right = size(spec_current_right, 2);
    
    CO_col_l = zeros(IM_left,  JM_left);
    CO_col_r = zeros(IM_right, JM_right);

    % bxheight [m] approximation only for now.
    bxheight_approx_l = zeros(IM_left, JM_left, LM_left, 12);
    bxheight_approx_r = zeros(IM_right, JM_right, LM_right, 12);
    for mo = 1:12
        bxheight_approx_l(:,:,:,mo) = bxheight_approx(file_left.T(:,:,:,mo), pedge_left);
        bxheight_approx_r(:,:,:,mo) = bxheight_approx(file_right.T(:,:,:,mo), pedge_right);
    end
    
    % 298K / 101325 Pa * 2.46e13 molec/cm3/vmr *
    %  * 100 cm/m
    phys_constant_factor = 298 / 101325 * 2.46e19 * 100;
    for L = 1:LM_left
        if fourth_idx > 0
            % molec/cm2
            CO_col_l = CO_col_l + spec_current_left(:,:,L,fourth_idx) .* (pedge_left(:,:,L) + pedge_left(:,:,L+1)) / 2 ./ file_left.T(:,:,L,fourth_idx) .* bxheight_approx_l(:,:,L,fourth_idx) * phys_constant_factor;
            %fprintf("debug: 25,25,%d, %f, %f, %f\n",L,spec_current_left(:,:,L,fourth_idx), file_left.T(25,25,L,fourth_idx), CO_col_l);
            %fprintf("debug: L=%d, spec=%e, bxheight_a=%e, sum=%e\n", L, spec_current_left(25,25,L,fourth_idx), bxheight_approx_l(25,25,L,fourth_idx), CO_col_l(25,25));
        elseif fourth_idx == -1
            for mo = 1:12
                CO_col_l = CO_col_l + spec_current_left(:,:,L,mo) .* (pedge_left(:,:,L) + pedge_left(:,:,L+1)) / 2 ./ file_left.T(:,:,L,mo) .* bxheight_approx_l(:,:,L,mo) * phys_constant_factor;
            end
            CO_col_l = CO_col_l / 12;
        end
    end

    for L = 1:LM_right
        if fourth_idx > 0
            % mol/cm2
            CO_col_r = CO_col_r + spec_current_right(:,:,L,fourth_idx) .* (pedge_right(:,:,L) + pedge_right(:,:,L+1)) / 2 ./ file_right.T(:,:,L,fourth_idx) .* bxheight_approx_r(:,:,L,fourth_idx) * phys_constant_factor;
        elseif fourth_idx == -1
            for mo = 1:12
                CO_col_r = CO_col_r + spec_current_right(:,:,L,mo) .* (pedge_right(:,:,L) + pedge_right(:,:,L+1)) / 2 ./ file_right.T(:,:,L,mo) .* bxheight_approx_r(:,:,L,mo) * phys_constant_factor;
            end
            CO_col_r = CO_col_r / 12;
        end
    end
    
    % sneakily auto fix cesm-circshift for CO_col_l and CO_col_r
    % handling of coordinates ...
    % cesm-circshift
    if max(lons_left) > 180.0
       CO_col_l = circshift(CO_col_l, ...
                     size(lons_left, 1)/2, ...
                     1); % data needs to be in (lat, lon) 

       lons_left_sf = lons_left - 180.0;

       fprintf("cesm-circshift: shifted data_left, lons_left\n");
    else
       lons_left_sf = lons_left;
    end

    if max(lons_right) > 180.0
       CO_col_r = circshift(CO_col_r, ...
                     size(lons_right, 1)/2, ...
                     1); % data needs to be in (lat, lon) 

       lons_right_sf = lons_right - 180.0;

       fprintf("cesm-circshift: shifted data_right, lons_right\n");
    else
       lons_right_sf = lons_right;
    end

    % also fix mopitt lats,lons to be 2D
    if size(mopitt_lons, 2) == 1
        lons_new = zeros(size(mopitt_lons, 1), size(mopitt_lats, 1));
        lats_new = zeros(size(mopitt_lons, 1), size(mopitt_lats, 1));
    
        for j = 1:size(mopitt_lats, 1)
            lons_new(:,j) = mopitt_lons(:,1);
        end
    
        for i = 1:size(mopitt_lons, 1)
            lats_new(i,:) = mopitt_lats(:,1);
        end
    
        mopitt_lons = lons_new;
        mopitt_lats = lats_new;
    end
    
    % plot column
    % LEFT
    axg = subplot_tight(2,3,1);

    hold on;
    m_proj('miller', ...
           'long',[-180.0 180.0], ...
           'lat', [-87.5 87.5]);

    m_pcolor(lons_left_sf, lats_left, CO_col_l);
    colorbar('southoutside');
    caxis(mopitt_co_range);

    m_plotbndry('/usr/local/MATLAB/R2021b/toolbox/boundary/world', ...
        'color', 'k', 'linewidth', 1.5)
    m_grid;
    colormap(axg, flip(othercolor('Spectral6', 12)));

    t = title(sprintf(left_name));
    set(gca, 'FontName', 'Arial');
    set(gca, 'FontSize', 16);
    set(findall(gcf,'type','text'), 'FontName', 'Arial');
    set(findall(gcf,'type','text'), 'FontSize', 16);

    % RIGHT
    axg = subplot_tight(2,3,2);

    hold on;
    m_proj('miller', ...
           'long',[-180.0 180.0], ...
           'lat', [-87.5 87.5]);

    m_pcolor(lons_right_sf, lats_right, CO_col_r);
    colorbar('southoutside');
    caxis(mopitt_co_range);

    m_plotbndry('/usr/local/MATLAB/R2021b/toolbox/boundary/world', ...
        'color', 'k', 'linewidth', 1.5)
    m_grid;
    colormap(axg, flip(othercolor('Spectral6', 12)));

    t = title(sprintf(right_name));
    set(gca, 'FontName', 'Arial');
    set(gca, 'FontSize', 16);
    set(findall(gcf,'type','text'), 'FontName', 'Arial');
    set(findall(gcf,'type','text'), 'FontSize', 16);
    
    
    % MOPITT CO observational data column [mol/cm2]
    axg = subplot_tight(2,3,3);

    hold on;
    m_proj('miller', ...
           'long',[-180.0 180.0], ...
           'lat', [-87.5 87.5]);

    m_pcolor(mopitt_lons, mopitt_lats, mopitt_co_column);
    colorbar('southoutside');
    caxis(mopitt_co_range);

    m_plotbndry('/usr/local/MATLAB/R2021b/toolbox/boundary/world', ...
        'color', 'k', 'linewidth', 1.5)
    m_grid;
    colormap(axg, flip(othercolor('Spectral6', 12)));

    t = title(sprintf("MOPITT CO observational data (Day) - %04d/%02d", mopitt_year, fourth_idx));
    set(gca, 'FontName', 'Arial');
    set(gca, 'FontSize', 16);
    set(findall(gcf,'type','text'), 'FontName', 'Arial');
    set(findall(gcf,'type','text'), 'FontSize', 16);
    
    data_l2m = griddata(lons_left_sf, lats_left, CO_col_l, ...
                        mopitt_lons, mopitt_lats);
    data_r2m = griddata(lons_right_sf, lats_right, CO_col_r, ...
                        mopitt_lons, mopitt_lats);
    delta_l2m = data_l2m - mopitt_co_column;
    delta_r2m = data_r2m - mopitt_co_column;
    delta_l2r = data_l2m - data_r2m;
    
    axg = subplot_tight(2,3,4);

    hold on;
    m_proj('miller', ...
           'long',[-180.0 180.0], ...
           'lat', [-87.5 87.5]);

    m_pcolor(mopitt_lons, mopitt_lats, delta_l2m);
    colorbar('southoutside');
    caxis(mopitt_co_delta_range);

    m_plotbndry('/usr/local/MATLAB/R2021b/toolbox/boundary/world', ...
        'color', 'k', 'linewidth', 1.5)
    m_grid;
    colormap(axg, flip(othercolor('RdBu11', 9)));

    t = title(sprintf("%s - MOPITT Col", left_name));
    set(gca, 'FontName', 'Arial');
    set(gca, 'FontSize', 16);
    set(findall(gcf,'type','text'), 'FontName', 'Arial');
    set(findall(gcf,'type','text'), 'FontSize', 16);
    
    
    axg = subplot_tight(2,3,5);

    hold on;
    m_proj('miller', ...
           'long',[-180.0 180.0], ...
           'lat', [-87.5 87.5]);

    m_pcolor(mopitt_lons, mopitt_lats, delta_r2m);
    colorbar('southoutside');
    caxis(mopitt_co_delta_range);


    m_plotbndry('/usr/local/MATLAB/R2021b/toolbox/boundary/world', ...
        'color', 'k', 'linewidth', 1.5)
    m_grid;
    colormap(axg, flip(othercolor('RdBu11', 9)));

    t = title(sprintf("%s - MOPITT Col", right_name));
    set(gca, 'FontName', 'Arial');
    set(gca, 'FontSize', 16);
    set(findall(gcf,'type','text'), 'FontName', 'Arial');
    set(findall(gcf,'type','text'), 'FontSize', 16);

    % Also do a diff of L-R while at it
    axg = subplot_tight(2,3,6);

    hold on;
    m_proj('miller', ...
           'long',[-180.0 180.0], ...
           'lat', [-87.5 87.5]);

    m_pcolor(mopitt_lons, mopitt_lats, delta_l2r);
    colorbar('southoutside');
    caxis(mopitt_co_delta_range);

    m_plotbndry('/usr/local/MATLAB/R2021b/toolbox/boundary/world', ...
        'color', 'k', 'linewidth', 1.5)
    m_grid;
    colormap(axg, flip(othercolor('RdBu11', 9)));

    t = title(sprintf("%s - %s (Diff L-R)", left_name, right_name));
    
    sgtitle(sprintf("CO column density [molec cm^{-2}] - date slice %d", fourth_idx));

    set(gcf, 'Renderer', 'painters', 'Position', [90 0 1600 1150]);
    set(gca, 'FontName', 'Arial');
    set(gca, 'FontSize', 16);
    set(findall(gcf,'type','text'), 'FontName', 'Arial');
    set(findall(gcf,'type','text'), 'FontSize', 16);
    
    fname = sprintf("%s_%02d_co_mopitt_column.png", figure_out_prefix, fourth_idx);
    saveas(gcf, fname);
    close(gcf);

    %-------------------------------------------------------------
    % also compare surface observations for CO against MOPITT.
    mopitt_co_range = [0, 300];
    mopitt_co_delta_range = [-150, 150];
    
    mopitt_co_sfc_ppb = transpose(mopitt_co_sfc_ppb_night); % huh? use day or night. fixme
    
    if fourth_idx > 0
        CO_sfc_l = file_left.CO(:,:,1,fourth_idx) * 1e9; % vmr to ppbv.
        CO_sfc_r = file_right.CO(:,:,1,fourth_idx) * 1e9;
    else
        CO_sfc_l = zeros(IM_left, JM_left);
        CO_sfc_r = zeros(IM_left, JM_left);
        for mo = 1:12
            CO_sfc_l = CO_sfc_l + file_left.CO(:,:,1,mo) * 1e9; % vmr to ppbv.
            CO_sfc_r = CO_sfc_l + file_right.CO(:,:,1,mo) * 1e9;
        end
    end
    
    % sneakily auto fix cesm-circshift for CO_sfc_l and CO_sfc_r
    % handling of coordinates ...
    % cesm-circshift
    if max(lons_left) > 180.0
       CO_sfc_l = circshift(CO_sfc_l, ...
                     size(lons_left, 1)/2, ...
                     1); % data needs to be in (lat, lon) 

       lons_left_sf = lons_left - 180.0;

       fprintf("cesm-circshift: shifted data_left, lons_left\n");
    else
       lons_left_sf = lons_left;
    end

    if max(lons_right) > 180.0
       CO_sfc_r = circshift(CO_sfc_r, ...
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

    m_pcolor(lons_left_sf, lats_left, CO_sfc_l);
    colorbar('southoutside');
    caxis(mopitt_co_range);

    m_plotbndry('/usr/local/MATLAB/R2021b/toolbox/boundary/world', ...
        'color', 'k', 'linewidth', 1.5)
    m_grid;
    colormap(axg, flip(othercolor('Spectral6', 12)));

    t = title(sprintf(left_name));
    set(gca, 'FontName', 'Arial');
    set(gca, 'FontSize', 16);
    set(findall(gcf,'type','text'), 'FontName', 'Arial');
    set(findall(gcf,'type','text'), 'FontSize', 16);

    % RIGHT
    axg = subplot_tight(2,3,2);

    hold on;
    m_proj('miller', ...
           'long',[-180.0 180.0], ...
           'lat', [-87.5 87.5]);

    m_pcolor(lons_right_sf, lats_right, CO_sfc_r);
    colorbar('southoutside');
    caxis(mopitt_co_range);

    m_plotbndry('/usr/local/MATLAB/R2021b/toolbox/boundary/world', ...
        'color', 'k', 'linewidth', 1.5)
    m_grid;
    colormap(axg, flip(othercolor('Spectral6', 12)));

    t = title(sprintf(right_name));
    set(gca, 'FontName', 'Arial');
    set(gca, 'FontSize', 16);
    set(findall(gcf,'type','text'), 'FontName', 'Arial');
    set(findall(gcf,'type','text'), 'FontSize', 16);
    
    
    % MOPITT CO observational data column [mol/cm2]
    axg = subplot_tight(2,3,3);

    hold on;
    m_proj('miller', ...
           'long',[-180.0 180.0], ...
           'lat', [-87.5 87.5]);

    m_pcolor(mopitt_lons, mopitt_lats, mopitt_co_sfc_ppb);
    colorbar('southoutside');
    caxis(mopitt_co_range);

    m_plotbndry('/usr/local/MATLAB/R2021b/toolbox/boundary/world', ...
        'color', 'k', 'linewidth', 1.5)
    m_grid;
    colormap(axg, flip(othercolor('Spectral6', 12)));

    t = title(sprintf("MOPITT CO retrieval (surface) - %04d/%02d", mopitt_year, fourth_idx));
    set(gca, 'FontName', 'Arial');
    set(gca, 'FontSize', 16);
    set(findall(gcf,'type','text'), 'FontName', 'Arial');
    set(findall(gcf,'type','text'), 'FontSize', 16);
    
    data_l2m = griddata(lons_left_sf, lats_left, CO_sfc_l, ...
                        mopitt_lons, mopitt_lats);
    data_r2m = griddata(lons_right_sf, lats_right, CO_sfc_r, ...
                        mopitt_lons, mopitt_lats);
    delta_l2m = data_l2m - mopitt_co_sfc_ppb;
    delta_r2m = data_r2m - mopitt_co_sfc_ppb;
    delta_l2r = data_l2m - data_r2m;
    
    axg = subplot_tight(2,3,4);

    hold on;
    m_proj('miller', ...
           'long',[-180.0 180.0], ...
           'lat', [-87.5 87.5]);

    m_pcolor(mopitt_lons, mopitt_lats, delta_l2m);
    colorbar('southoutside');
    caxis(mopitt_co_delta_range);

    m_plotbndry('/usr/local/MATLAB/R2021b/toolbox/boundary/world', ...
        'color', 'k', 'linewidth', 1.5)
    m_grid;
    colormap(axg, flip(othercolor('RdBu11', 9)));

    t = title(sprintf("%s - MOPITT Surface VMR", left_name));
    set(gca, 'FontName', 'Arial');
    set(gca, 'FontSize', 16);
    set(findall(gcf,'type','text'), 'FontName', 'Arial');
    set(findall(gcf,'type','text'), 'FontSize', 16);
    
    
    axg = subplot_tight(2,3,5);

    hold on;
    m_proj('miller', ...
           'long',[-180.0 180.0], ...
           'lat', [-87.5 87.5]);

    m_pcolor(mopitt_lons, mopitt_lats, delta_r2m);
    colorbar('southoutside');
    caxis(mopitt_co_delta_range);


    m_plotbndry('/usr/local/MATLAB/R2021b/toolbox/boundary/world', ...
        'color', 'k', 'linewidth', 1.5)
    m_grid;
    colormap(axg, flip(othercolor('RdBu11', 9)));

    t = title(sprintf("%s - MOPITT Surface VMR", right_name));
    set(gca, 'FontName', 'Arial');
    set(gca, 'FontSize', 16);
    set(findall(gcf,'type','text'), 'FontName', 'Arial');
    set(findall(gcf,'type','text'), 'FontSize', 16);

    % Also do a diff of L-R while at it
    axg = subplot_tight(2,3,6);

    hold on;
    m_proj('miller', ...
           'long',[-180.0 180.0], ...
           'lat', [-87.5 87.5]);

    m_pcolor(mopitt_lons, mopitt_lats, delta_l2r);
    colorbar('southoutside');
    caxis(mopitt_co_delta_range);

    m_plotbndry('/usr/local/MATLAB/R2021b/toolbox/boundary/world', ...
        'color', 'k', 'linewidth', 1.5)
    m_grid;
    colormap(axg, flip(othercolor('RdBu11', 9)));

    t = title(sprintf("%s - %s (Diff L-R)", left_name, right_name));
    
    sgtitle(sprintf("CO surface conc [ppbv] - date slice %d", fourth_idx));

    set(gcf, 'Renderer', 'painters', 'Position', [90 0 1600 1150]);
    set(gca, 'FontName', 'Arial');
    set(gca, 'FontSize', 16);
    set(findall(gcf,'type','text'), 'FontName', 'Arial');
    set(findall(gcf,'type','text'), 'FontSize', 16);
    
    fname = sprintf("%s_%02d_co_mopitt_sfc.png", figure_out_prefix, fourth_idx);
    saveas(gcf, fname);
    close(gcf);
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
% for the entire grid. [Pa]
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

% bxheight_approx: gives approximation of box height [m] for each
% grid box. this is an approximation because we do not use
% moist air in this calculation, instead assuming q = 0.
% virtual temperature should be used instead, but Q/H2O is not always
% available.
function bxheight_approx = bxheight_approx(T, PEDGE)
    IM = size(PEDGE, 1);
    JM = size(PEDGE, 2);
    LM = size(PEDGE, 3) - 1;

    bxheight_approx = zeros(IM, JM, LM);
    for L = 1:LM
        bxheight_approx(:,:,L) = (287.0 / 9.80665) * T(:,:,L) .* log(PEDGE(:,:,L) ./ PEDGE(:,:,L+1));
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

% Changelog:
% 2022.09.15 - Quality of life improvements. Allow reference_model to read
% already computed data for comparison.
%
