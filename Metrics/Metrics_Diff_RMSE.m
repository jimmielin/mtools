%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%            _____            ______        %
% _______ _____  /_______________  /_______ %
% __  __ `__ \  __/  __ \  __ \_  /__  ___/ %
% _  / / / / / /_ / /_/ / /_/ /  / _(__  )  %
% /_/ /_/ /_/\__/ \____/\____//_/  /____/   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       "mtools" Research Toolkit           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Metrics_Diff_RMSE.m
%
% Computes RMSE differences for two model runs
% grid-box-by-grid-box (same grid required)
%
% Now also:
% - Computes RRMS differences
%   > all species, all domain
%   > by species
%   > by vertical level
% - Creates RRMS 
%   > histograms,
%   > by-level plots, 
%   > vertically integrated plots by species,
%   > by-species totals, 
%   > totals     
% - Allows for
%   > geographical subsetting
%   > defined timestep (run from external loop)
%   >> if defined timestep, log histogram peaks if histogram
% - with Threshold to avoid small-difference-small-value problem
%   > based on absolute molec/cm3 threshold
%   > based on relative threshold by layer (i.e., 10% of mean)
% - can be driven by external MATLAB prompt through external_timestr
%   > ... and save 'all', 'perspc', and 'histogram' peaks into external_save_.. vars
%   > ... and save median, 5/25/75/95th percentile of all spc median RRMS to external vars
% - merges with regard to isorropia balancing, two families of species:
%   NH3T = NH3 + NH4, and HNO3T = HNO3 + NIT
%
% Reads raw netCDF files for convenience
%
% (c) 2019-2022 Haipeng Lin <jimmie.lin@gmail.com>
% Cambridge, MA | Washington, DC | Munich | Paris | Lisbon
%
% Version: 2022.04.20
% Started: 2021.11.16
% See changelog at end of script%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% specify data on left and right
%
% netCDF file names (standard model output)
%
% RMSE percentages (RRMS) are calculated based on filename_left.
casename = 'expX_5.1-4';
base_casename = 'exp1_baseline';

if exist('external_timestr', 'var')
    timestr = external_timestr;
else
    timestr  = '20150212_0000';
end

filename_left = sprintf('%s/GEOSChem.SpeciesConcMND.%sz.nc4', base_casename, timestr);
filename_right = sprintf('%s/GEOSChem.SpeciesConcMND.%sz.nc4', casename, timestr);

%% below file-independent code (separated for optimizing section run)
% specify variables to calculate RMSE/RRMS for.
%
% if mode == 'histogram', specify vars_in_3D as 9 species for a 3x3 grid

vars_in_2D = [];
%vars_in_3D = ["SpeciesConcMND_O3", "SpeciesConcMND_NO", "SpeciesConcMND_NO2", "SpeciesConcMND_CH2O", "SpeciesConcMND_OH", "SpeciesConcMND_CO", "SpeciesConcMND_ISOP", "SpeciesConcMND_SO2", "SpeciesConcMND_SO4", "SpeciesConcMND_NIT", ...
%    "SpeciesConcMND_HO2"];

% some families:
% HOx:
% ISORROPIA:
% vars_in_3D = ["SpeciesConcMND_HNO3", "SpeciesConcMND_NH3", "SpeciesConcMND_NH4", "SpeciesConcMND_NIT", "SpeciesConcMND_SALA", "SpeciesConcMND_SO4", "SpeciesConcMND_SALACL", "SpeciesConcMND_HCl", "SpeciesConcMND_SALC", "SpeciesConcMND_SALCCL", "SpeciesConcMND_NITs", "SpeciesConcMND_SO4s", "SpeciesConcMND_SALAAL", "SpeciesConcMND_SALCAL"];

% Halogens:
% vars_in_3D = ["SpeciesConcMND_Cl", "SpeciesConcMND_ClO", "SpeciesConcMND_HCl", "SpeciesConcMND_ICl", "SpeciesConcMND_BrCl", "SpeciesConcMND_ClNO3", "SpeciesConcMND_Cl2", "SpeciesConcMND_HOCl", "SpeciesConcMND_ClNO2", "SpeciesConcMND_Br", "SpeciesConcMND_BrO", "SpeciesConcMND_HBr", "SpeciesConcMND_IBr", "SpeciesConcMND_BrNO3", "SpeciesConcMND_Br2", "SpeciesConcMND_HOBr", "SpeciesConcMND_BrSALA", "SpeciesConcMND_BrSALC", "SpeciesConcMND_BrNO2", "SpeciesConcMND_CH2IBr"];

% main and species of interest
% vars_in_3D = ["SpeciesConcMND_O3", "SpeciesConcMND_NO", "SpeciesConcMND_NO2", "SpeciesConcMND_CH2O", "SpeciesConcMND_OH", "SpeciesConcMND_CO", "SpeciesConcMND_SO2", "SpeciesConcMND_SALA", "SpeciesConcMND_SO4", "SpeciesConcMND_NITs", "SpeciesConcMND_SO4s", "SpeciesConcMND_Cl", "SpeciesConcMND_ClO", "SpeciesConcMND_ICl", "SpeciesConcMND_BrCl", "SpeciesConcMND_ClNO3", "SpeciesConcMND_Cl2", "SpeciesConcMND_HOCl", "SpeciesConcMND_ClNO2", "SpeciesConcMND_Br", "SpeciesConcMND_BrO", "SpeciesConcMND_HBr", "SpeciesConcMND_IBr", "SpeciesConcMND_BrNO3", "SpeciesConcMND_Br2", "SpeciesConcMND_HOBr", "SpeciesConcMND_BrNO2", "SpeciesConcMND_CH2IBr", "SpeciesConcMND_NH3T", "SpeciesConcMND_HNO3T"];

% vars_in_3D = ["SpeciesConcMND_I2", "SpeciesConcMND_INO"];

% - 9 spc list for plotting histograms -
% Halogens and main:
% vars_in_3D = ["SpeciesConcMND_Br", "SpeciesConcMND_BrO", "SpeciesConcMND_HBr", ...
%               "SpeciesConcMND_Cl", "SpeciesConcMND_ClO", "SpeciesConcMND_HCl", ...
%               "SpeciesConcMND_O3", "SpeciesConcMND_NO2", "SpeciesConcMND_CO"];

% ISORROPIA spc:
% vars_in_3D = ["SpeciesConcMND_NH3", "SpeciesConcMND_NH4", "SpeciesConcMND_NIT", ...
%               "SpeciesConcMND_HNO3", "SpeciesConcMND_SO2", "SpeciesConcMND_SO4", ...
%               "SpeciesConcMND_SALA", "SpeciesConcMND_SALC", "SpeciesConcMND_NITs"];

% debug acids:
%vars_in_3D = ["SpeciesConcMND_NO2", "SpeciesConcMND_SO4s", "SpeciesConcMND_SO4", ...
%              "SpeciesConcMND_HNO3", "SpeciesConcMND_HCl", "SpeciesConcMND_SALA", ...
%              "SpeciesConcMND_SALC", "SpeciesConcMND_BrSALA", "SpeciesConcMND_BrSALC"];

vars_in_3D = ["SpeciesConcMND_I", "SpeciesConcMND_I2O3", "SpeciesConcMND_HI", ...
              "SpeciesConcMND_I2O2", "SpeciesConcMND_OIO", "SpeciesConcMND_ICl", ...
              "SpeciesConcMND_IONO2", "SpeciesConcMND_NITs", "SpeciesConcMND_SO4s"]; % testing only

% mode: all|all-lev|perspc|perspc-autoall|split|histogram|map
%
% all:     compute RMSE over the entire vertical for ALL 3-D SpeciesConcMND_
% all-lev: compute RMSE over each level individually for ALL 3-D SpeciesConcMND_
% perspc:  compute RMSE over the entire vertical for some vars
%  perspc-autoall: all, but no need to fill vars_in_3D yourself
% perspc2: compute RMSE over the entire specified vertical for some vars.
%          uses a masking method from Shen et al., 2022, where 99% of the mass will be accounted for
%          and the remaining grid boxes dropped.
%  perspc2-autoall: same
%
% split:   compute RMSE over each level individually plotting end result, for some vars
% histogram: RRMS over every grid box and plot histogram
% map:     make a m_map (boundary package required) for RRMS, vertically integrated
mode = "perspc2-autoall";

% subset for perspc2: sfc|500|trop|strat|pbl|freetrop|all
% now uses climatological tropopause heights
perspc2_subset = "all";

% tropopause level climatology data path format string (%02d = month)
% has to supply 1-year monthly tropopause levels for use in perspc2_subset
troplev_climo_data = "/n/holyscratch01/jacob_lab/hplin/gc_dev_202203/tropopause_climatology/OutputDir/GEOSChem.StateMet.2014%02d01_0000z.nc4";
timestr_month = str2num(timestr(5:6));
troplev_climo_data_file = sprintf(troplev_climo_data, timestr_month);

% debug only for perspc - prints out per layer statistics
perspc_debug = 0;

% perspc_flag: custom features
%
% 0: none
% 1: plot RRMS % against mean concentration (requires extensive extra compute)
% -- calculate RRMS % and split by strat-trop (useful for externally forced) --
% 5: 500 hPa only
% 6: surface only
% 7: ..... strat only
% 8: ..... trop only
perspc_flag = 0;

% use_median_not_mean: 0|1
% over the I,J,L dimensions use median as a metric and not mean of the errors
% (formerly perspc_flag 9)
use_median_not_mean = 0;

% histogram_scale: log | regular (if mode == 'histogram')
histogram_scale = "log";

% mask_method: relative|absolute
%
% if absolute, any value_left < mask_threshold will be ignored
mask_method = "relative";

% mask_threshold (if mask_method == 'absolute')
% set to 1e6 molec/cm3 per Santillana et al., 2010; Shen et al., 2020
%
% since data is in v/v dry, convert ppv to molec/cm3 roughly using
% 1 atm ~ 2.46e19 molec/cm3
%
%
% used to avoid division by zero
mask_threshold = 10; % 1e4

% if relative,
%   if rel_ratio > 0,
%   any value_left < rel_ratio * vertical layer mean will be ignored
%   rel_ratio usually chosen 0.1 (ignore lower than 10% * avg) or 0.01
%            , AND any value_left < mask_threshold will be ignored
%
%   if rel_ratio < 0, then instead filters out -rel_ratio (th) percentile
%   e.g., -1 will remove grid boxes corresp. lowest 1% percentile of this level
%   using left_model as a reference
mask_threshold_rel_ratio = 0.01;

% IM, JM subsetting
% do any subsetting of data? specify range if so.
%
% ignored for world plots obviously
%
% warning: no effort is made to check IM_, JM_range validity
horiz_subset = false;
IM_range = 13:17;
JM_range = 25:29;

% world_boundary_file: '/usr/local/MATLAB/R2021b/toolbox/boundary/world'
% get this from the boundary toolbox.
world_boundary_file = '/n/home04/hplin/boundary/world';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% NO USER CONFIGURABLE CODE BELOW %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf("Comparing:\nLeft_Model: %s\nRight_Model: %s\n", filename_left, filename_right);
fprintf("Mode: %s\n", mode);

if strcmp(mode, "perspc2")
    fprintf("=> Perspc2 mode: %s\n", perspc2_subset);
end

% if all species, replace vars_in_3D with all the SpeciesConcMND_
%
isautoallflag = false;
if strcmp(mode, "all") || strcmp(mode, "perspc-autoall") || strcmp(mode, "perspc2-autoall")
    vars_in_3D = ["dummy"]; % because otherwise we cannot define an all-string list...
    vars_raw_nc = {ncinfo(filename_left).Variables.Name};
    for i = 1:length(vars_raw_nc)
        if contains(char(vars_raw_nc(i)), "SpeciesConcMND_")
            % fprintf("check: %s\n", char(vars_raw_nc(i)))

            % Skip HNO3, NH3, NH4, NIT and merge
            % if contains(char(vars_raw_nc(i)), "SpeciesConcMND_HNO3") || contains(char(vars_raw_nc(i)), "SpeciesConcMND_NIT") || contains(char(vars_raw_nc(i)), "SpeciesConcMND_NH3") || contains(char(vars_raw_nc(i)), "SpeciesConcMND_NH4")
            %     if ~contains(char(vars_raw_nc(i)), "NITs")
            %         continue
            %     end
            % end
            %
            % does not skip anymore so they can be accessed

            % Blacklist species which have no valid concentrations (all NaN/-999)
            % these are purely "dummy species" as defined by the species database
            if contains(char(vars_raw_nc(i)), "LNRO2H") || contains(char(vars_raw_nc(i)), "LNRO2N") || contains(char(vars_raw_nc(i)), "NAP") || contains(char(vars_raw_nc(i)), "NRO2") || contains(char(vars_raw_nc(i)), "SO4H1") || contains(char(vars_raw_nc(i)), "SO4H2") || contains(char(vars_raw_nc(i)), "SO4H3") || contains(char(vars_raw_nc(i)), "SO4H4") || contains(char(vars_raw_nc(i)), "LISOPNO3") || contains(char(vars_raw_nc(i)), "LISOPOH") || contains(char(vars_raw_nc(i)), "LNRO2H") || contains(char(vars_raw_nc(i)), "LNRO2N") || contains(char(vars_raw_nc(i)), "LOx") || contains(char(vars_raw_nc(i)), "LTRO2H") || contains(char(vars_raw_nc(i)), "LTRO2N") || contains(char(vars_raw_nc(i)), "LXRO2H") || contains(char(vars_raw_nc(i)), "LXRO2N") || contains(char(vars_raw_nc(i)), "PCO") || contains(char(vars_raw_nc(i)), "PH2O2") || contains(char(vars_raw_nc(i)), "POx") || contains(char(vars_raw_nc(i)), "LBRO2H") || contains(char(vars_raw_nc(i)), "LBRO2N") || contains(char(vars_raw_nc(i)), "LCH4") || contains(char(vars_raw_nc(i)), "LCO") || contains(char(vars_raw_nc(i)), "PSO4")
                continue
            end

            vars_in_3D(end+1) = sprintf("%s", char(vars_raw_nc(i)));
        end
    end
    vars_in_3D(end+1) = "SpeciesConcMND_NH3T";
    vars_in_3D(end+1) = "SpeciesConcMND_HNO3T";
    vars_in_3D(end+1) = "SpeciesConcMND_HNO3Ts";
    vars_in_3D(end+1) = "SpeciesConcMND_SO4Ts";
    vars_in_3D(1) = []; % remove the dummy element

    if strcmp(mode, "perspc-autoall")
        mode = "perspc";
        isautoallflag = true;
    end

    if strcmp(mode, "perspc2-autoall")
        mode = "perspc2";
        isautoallflag = true;
    end
end

% sort the species list in order
if ~strcmp(mode, "histogram")
    vars_in_3D = sort(vars_in_3D);
end


% --------
% read the coords
% fixme: assuming same grid
lons = ncread(filename_left, 'lon');
lats = ncread(filename_left, 'lat');
levs = ncread(filename_left, 'lev');

IM   = numel(lons);
JM   = numel(lats);
LM   = numel(levs);

% if horiz_subset is set to true, then IM, JM total dims will change
if horiz_subset
    IM = length(IM_range);
    JM = length(JM_range);

    ips = IM_range(1);
    ipe = IM_range(end);
    jps = JM_range(1);
    jpe = JM_range(end);

    fprintf("Using horizontal subsetting. ips, ipe, jps, jpe: %d, %d, %d, %d\n", ips, ipe, jps, jpe);
else
    ips = 1;
    ipe = IM;
    jps = 1;
    jpe = JM;
end

if use_median_not_mean
    fprintf("Using median for geographical metrics, not mean.\n");
end

% check if we need to flip the vertical. this is so that TOA is always the
% highest indexed level, by convention in mtools
if levs(1) < levs(2)
    levs = flip(levs);
    do_flip_vert = true;
    fprintf("Notice: levels were found as TOA = level 1 in input file. By convention in this set of tools, the vertical has been flipped so that coordinates face up.\n");
else
    do_flip_vert = false;
end

% check if we need to make 1-D rectilinear coordinates into a 2-D one
% for mapping...
if size(lons, 2) == 1
    lons_2d = zeros(size(lons, 1), size(lats, 1));
    lats_2d = zeros(size(lons, 1), size(lats, 1));

    for j = 1:size(lats, 1)
        lons_2d(:,j) = lons(:,1);
    end

    for i = 1:size(lons, 1)
        lats_2d(i,:) = lats(:,1);
    end
end

% if perspc, store rrms per species so we can sort later
%
% also compute quantities for tropopause and PBL height levels
if strcmp(mode, "perspc") || strcmp(mode, "perspc2")
    perspc_rrms = zeros(length(vars_in_3D), 1);
    perspc_rmse = zeros(length(vars_in_3D), 1);
    perspc_rrms_byName = containers.Map();

    %---------------------------------------------
    % eventually read the tropopause height mask.
    %---------------------------------------------
    trop_level = ncread(troplev_climo_data_file, "Met_TropLev");
    trop_level = round(trop_level); % this is monthly data...

    pbl_level  = ncread(troplev_climo_data_file, "Met_PBLTOPL");
    pbl_level  = round(pbl_level);

    read_layer_ks = ones(IM, JM);
    read_layer_ke = zeros(IM, JM);
    % subset for perspc2: sfc|500|trop|strat
    if strcmp(perspc2_subset, "sfc")
        read_layer_ke = ones(IM, JM);
    elseif strcmp(perspc2_subset, "500")
        read_layer_ks = ones(IM, JM) * 23;
        read_layer_ke = ones(IM, JM) * 23;
    elseif strcmp(perspc2_subset, "pbl")
        read_layer_ke = pbl_level;
    elseif strcmp(perspc2_subset, "freetrop")
        read_layer_ks = pbl_level + 1;
        read_layer_ke = trop_level;
    elseif strcmp(perspc2_subset, "trop")
        % read_layer_ke = ones(IM, JM) * 35; % to read trop height.
        read_layer_ke = trop_level;
    elseif strcmp(perspc2_subset, "strat")
        % read_layer_ks = ones(IM, JM) * 36;
        read_layer_ks = trop_level + 1;
        read_layer_ke = ones(IM, JM) * LM;
    elseif strcmp(perspc2_subset, "all")
        read_layer_ke = ones(IM, JM) * LM;
    end

    % air volume data
    airvol_file = "/n/holyscratch01/jacob_lab/hplin/gc_dev_202203/tropopause_climatology/OutputDir2/GEOSChem.StateMet.20140101_0000z.nc4";
    airvol_3D = ncread(airvol_file, "Met_AIRVOL");
end

% if mode is in histogram, allow for "listening" for species to note their
% peak values.
%
% this is one of these "easter egg" hacky features ...
if strcmp(mode, "histogram") && exist('external_timestr', 'var')
    fprintf("external_timestr: a dictionary external_log_histogram_peak_trop|strat will be available to log histogram peaks...")
    if ~exist('external_log_histogram_peak_trop', 'var')
        fprintf(" (creating now)\n");
        external_log_histogram_peak_trop  = containers.Map();
        external_log_histogram_peak_strat = containers.Map();
        for i = 1:length(vars_in_3D)
            rawSpcName = strrep(vars_in_3D(i), "SpeciesConcMND_", "");
            external_log_histogram_peak_trop(rawSpcName) = [];
            external_log_histogram_peak_strat(rawSpcName) = [];
        end
    else
        fprintf(" will APPEND.\n");
    end
end

% if mode is all or perspc, also save external_log_rrms_all or
% external_log_rrms_perspc(spc)
if strcmp(mode, "all") && exist('external_timestr', 'var')
    fprintf("external_timestr: a dictionary external_log_rrms_all will be available to log total RRMS %...")
    if ~exist('external_log_rrms_all', 'var')
        fprintf(" (creating now)\n");
        external_log_rrms_all  = [];
    else
        fprintf(" will APPEND.\n");
    end
end

if (strcmp(mode, "perspc") || strcmp(mode, "perspc2")) && exist('external_timestr', 'var')
    fprintf("external_timestr: a dictionary external_log_rrms|rmse_perspc will be available to log per-species RRMS %...")
    if ~exist('external_log_rrms_perspc', 'var')
        fprintf(" (creating now)\n");
        external_log_rrms_perspc = containers.Map();
        external_log_rmse_perspc = containers.Map();
        external_log_val_left_mean_perspc = containers.Map();
        external_log_val_right_mean_perspc = containers.Map();

        if isautoallflag == true
            external_log_rrms_allspc_mean   = [];
            external_log_rrms_allspc_median = [];
            external_log_rrms_allspc_prc05  = [];
            external_log_rrms_allspc_prc25  = [];
            external_log_rrms_allspc_prc75  = [];
            external_log_rrms_allspc_prc95  = [];
        end

        for i = 1:length(vars_in_3D)
            rawSpcName = strrep(vars_in_3D(i), "SpeciesConcMND_", "");
            external_log_rrms_perspc(rawSpcName) = [];
            external_log_rmse_perspc(rawSpcName) = [];
            external_log_val_left_mean_perspc(rawSpcName) = [];
            external_log_val_right_mean_perspc(rawSpcName) = [];
        end
    else
        fprintf(" will APPEND.\n");
    end
end

% FIXME: using lev as rough estimate and not hy(a|b)(i|m) as we should
% but this is a decent approximation. also, vertical should match
% for both files anyway so it is just the labeling of prs which is not
% really accurate

% --------
% read variable one by one into buffer. note this reads left and right
% together, which means that if this is modified for diff grids in the
% future may need a loop split
for i = 1:length(vars_in_3D)
    spiename = vars_in_3D(i);
    if strcmp(spiename, "dummy")
        continue
    end
    rawSpcName = strrep(spiename, "SpeciesConcMND_", "");

    % we only support one subslice at this time, so always read the first
    % sub-slice of the file (TODO 12/6/21)

    % check if not special merged case...
    if strcmp(rawSpcName, "NH3T")
        tmp_left  = ncread(filename_left,  "SpeciesConcMND_NH3") + ncread(filename_left,  "SpeciesConcMND_NH4");
        tmp_right = ncread(filename_right, "SpeciesConcMND_NH3") + ncread(filename_right, "SpeciesConcMND_NH4");
    elseif strcmp(rawSpcName, "HNO3T")
        tmp_left  = ncread(filename_left,  "SpeciesConcMND_HNO3") + ncread(filename_left,  "SpeciesConcMND_NIT");
        tmp_right = ncread(filename_right, "SpeciesConcMND_HNO3") + ncread(filename_right, "SpeciesConcMND_NIT");
    elseif strcmp(rawSpcName, "HNO3Ts")
        tmp_left  = ncread(filename_left,  "SpeciesConcMND_HNO3") + ncread(filename_left,  "SpeciesConcMND_NIT") + ncread(filename_left,  "SpeciesConcMND_NITs");
        tmp_right = ncread(filename_right, "SpeciesConcMND_HNO3") + ncread(filename_right, "SpeciesConcMND_NIT") + ncread(filename_right, "SpeciesConcMND_NITs");
    elseif strcmp(rawSpcName, "SO4Ts")
        tmp_left  = ncread(filename_left,  "SpeciesConcMND_SO4") + ncread(filename_left,  "SpeciesConcMND_SO4s");
        tmp_right = ncread(filename_right, "SpeciesConcMND_SO4") + ncread(filename_right, "SpeciesConcMND_SO4s");
    else
        tmp_left  = ncread(filename_left,  spiename);
        tmp_right = ncread(filename_right, spiename);
    end
    superdim_left = size(size(tmp_left));

    if superdim_left(2) == 4
        tmp_left = tmp_left(:,:,:,1);
        fprintf("Warning: multiple time slices were found within tmp_left/right; only first time slice is used. Multiple slices in one file are presently unsupported\n");

        tmp_right = tmp_right(:,:,:,1);
    end

    if strcmp(mode, "perspc")
        if strcmp(mask_method, "relative")
            % load into chunks ...
            unified_chunk_size = IM*JM;

            % check if chunk size is custom. if so, OVERRIDE LM ranges...
            % strat only.
            if perspc_flag == 7
                tmp_left_1d  = zeros(unified_chunk_size*(LM-35+1), 1);
                tmp_right_1d = zeros(unified_chunk_size*(LM-35+1), 1);
                L_range = 35:LM;
                L_idx_relative = 35-1;
            elseif perspc_flag == 8 % trop only
                tmp_left_1d  = zeros(unified_chunk_size*35, 1);
                tmp_right_1d = zeros(unified_chunk_size*35, 1);
                L_range = 1:35;
                L_idx_relative = 1-1;
            elseif perspc_flag == 6 % surface only
                tmp_left_1d  = zeros(unified_chunk_size*1, 1);
                tmp_right_1d = zeros(unified_chunk_size*1, 1);
                L_range = 1:1;
                L_idx_relative = 1-1;
            elseif perspc_flag == 5 % 500 hPa only
                tmp_left_1d  = zeros(unified_chunk_size*1, 1);
                tmp_right_1d = zeros(unified_chunk_size*1, 1);
                L_range = 23:23;
                L_idx_relative = 23-1;
            else
                tmp_left_1d  = zeros(unified_chunk_size*LM, 1);
                tmp_right_1d = zeros(unified_chunk_size*LM, 1);
                L_range = 1:LM;
                L_idx_relative = 1-1;
            end

            for L = L_range
                tmp_left_2d = tmp_left(ips:ipe,jps:jpe,L);

                if mask_threshold_rel_ratio >= 0
                    l_avg_thres = mask_threshold_rel_ratio * mean(tmp_left_2d(:));
                else
                    % remove lowermost -ratio percentile ... such a badly named feature
                    l_avg_thres = prctile(tmp_left_2d(:), -mask_threshold_rel_ratio);
                end

                if l_avg_thres < 1e-6
                    % average threshold is too small - this whole layer is zero!
                    l_avg_thres = 1e-6;
                end

                if perspc_debug
                    perspc_debug_mean = mean(tmp_left_2d(:));
                    perspc_debug_median = median(tmp_left_2d(:));
                    perspc_debug_max = max(tmp_left_2d(:));
                    % do a trial compute (resource intensive)
                    tmp_left_2d(tmp_left_2d <= l_avg_thres) = -999;
                    tmp_right_2d = tmp_right(ips:ipe,jps:jpe,L);
                    tmp_right_2d(tmp_left_2d <= l_avg_thres) = -999;
                    perspc_debug_left_1d = tmp_left_2d(:);
                    perspc_debug_right_1d = tmp_right_2d(:);
                    perspc_debug_mask_1d = perspc_debug_left_1d > mask_threshold;
                    perspc_debug_rmse = sqrt(mean(((perspc_debug_right_1d(perspc_debug_mask_1d) - perspc_debug_left_1d(perspc_debug_mask_1d))) .^ 2));
                    perspc_debug_rrms = sqrt(mean(((perspc_debug_right_1d(perspc_debug_mask_1d) - perspc_debug_left_1d(perspc_debug_mask_1d))./perspc_debug_left_1d(perspc_debug_mask_1d)) .^ 2)) * 100;

                    fprintf("%s @ L = %d - mean %.4e, median %.4e, max %.4e, thres %.4e -- %.4f%s (rmse %.4f)\n", rawSpcName, L, perspc_debug_mean, perspc_debug_median, perspc_debug_max, l_avg_thres, perspc_debug_rrms, '%', perspc_debug_rmse);
                end

                % because a mask needs to be applied later, set -999 to use the same mask_threshold pass
                tmp_left_2d(tmp_left_2d <= l_avg_thres) = -999;

                tmp_right_2d = tmp_right(ips:ipe,jps:jpe,L);
                tmp_right_2d(tmp_left_2d <= l_avg_thres) = -999;

                % reshape into 1-D array...
                tmp_left_1d(1+unified_chunk_size*(L-1-L_idx_relative):unified_chunk_size*(L-L_idx_relative)) = reshape(tmp_left_2d, 1, unified_chunk_size);
                tmp_right_1d(1+unified_chunk_size*(L-1-L_idx_relative):unified_chunk_size*(L-L_idx_relative)) = reshape(tmp_right_2d, 1, unified_chunk_size);
            end

            mask = tmp_left_1d > mask_threshold;
            if isempty(tmp_left_1d(mask))
                fprintf("Warning: %s has no valid elements after threshold; before threshold min = %e\n", spiename, min(tmp_left_1d(:)));
                rmse = 0; rrms = 0;
            end
        elseif strcmp(mask_method, "absolute")
            mask = tmp_left_1d > mask_threshold;
            tmp_left_1d = tmp_left(:);
            tmp_right_1d = tmp_right(:);
        end

        if use_median_not_mean == 1
            rmse = sqrt(median((tmp_right_1d(mask) - tmp_left_1d(mask)) .^ 2));
            rrms = sqrt(median(((tmp_right_1d(mask) - tmp_left_1d(mask))./tmp_left_1d(mask)) .^ 2)) * 100;
        else
            rmse = sqrt(mean((tmp_right_1d(mask) - tmp_left_1d(mask)) .^ 2));
            rrms = sqrt(mean(((tmp_right_1d(mask) - tmp_left_1d(mask))./tmp_left_1d(mask)) .^ 2)) * 100;
        end

        perspc_rrms(i) = rrms;
        perspc_rmse(i) = rmse;
        perspc_rrms_byName(rawSpcName) = rrms;

        if exist('external_log_rrms_perspc', 'var')
            external_log_rrms_perspc(rawSpcName)  = [external_log_rrms_perspc(rawSpcName), rrms];
            external_log_rmse_perspc(rawSpcName)  = [external_log_rmse_perspc(rawSpcName), rmse];
            external_log_val_left_mean_perspc(rawSpcName) = [external_log_val_left_mean_perspc(rawSpcName), mean(tmp_left_1d(mask))];
            external_log_val_right_mean_perspc(rawSpcName) = [external_log_val_right_mean_perspc(rawSpcName), mean(tmp_right_1d(mask))];
        end

        % custom flag 1: plot individual spc rrms - concentration scatter plot
        if perspc_flag == 1
            % make temp 1d sizes
            mask_size = length(tmp_left_1d(mask));
            rrs = sqrt((tmp_right_1d(mask) - tmp_left_1d(mask) ./ tmp_left_1d(mask)) .^ 2) * 100;
            rse = sqrt((tmp_right_1d(mask) - tmp_left_1d(mask)) .^ 2);

            % plot rrs against tmp_left_1d in scatter plot ...
            scatter(tmp_left_1d(mask), rrs, 4, 'filled');
            xlabel("left_model value", "FontSize", 14, 'Interpreter', 'none');
            ylabel("relative error in grid box [%]", "FontSize", 14, 'Interpreter', 'none');
            title(sprintf("(rrms) Species = %s", rawSpcName), "FontSize", 14);
            subtitle(sprintf("%s vs %s", filename_left, filename_right), 'Interpreter', 'none', 'FontSize', 10);

            set(gcf, 'Renderer', 'painters', 'Position', [90 0 800 640]);

            % just save
            fname = sprintf("RRS_1_bySpc_%s_%s.png", timestr, rawSpcName);
            saveas(gcf, sprintf("out/%s", fname));
            close(gcf);

            % plot rse against tmp_left_1d in scatter plot ...
            scatter(tmp_left_1d(mask), rse, 4, 'filled');
            xlabel("left_model value", "FontSize", 14, 'Interpreter', 'none');
            ylabel("absolute error in grid box", "FontSize", 14, 'Interpreter', 'none');
            title(sprintf("(rmse) Species = %s", rawSpcName), "FontSize", 14);
            subtitle(sprintf("%s vs %s", filename_left, filename_right), 'Interpreter', 'none', 'FontSize', 10);

            set(gcf, 'Renderer', 'painters', 'Position', [90 0 800 640]);

            % just save
            fname = sprintf("RSE_1_bySpc_%s_%s.png", timestr, rawSpcName);
            saveas(gcf, sprintf("out/%s", fname));
            close(gcf);
        end

        fprintf("=> %s RMSE = %.4e, RRMS = %.4f%s (min=%.2e, avg=%.2e)\n", rawSpcName, rmse, rrms, '%', 0.0, 0.0);
    % --------------------- PERSPC2 (Shen et al., 2022 method) ---------------------
    elseif strcmp(mode, "perspc2")
        % the _ke and _ks indices are read outside of the loop for efficiency

        %---------------------------------------------
        % extract the data into a 1-dimensional chunk.
        %---------------------------------------------
        chunk_diff = read_layer_ke - read_layer_ks + 1; % this is a per-element sum.
        total_chunk_size = sum(chunk_diff, 'all');

        % allocate the memory.
        tmp_left_1d = zeros(total_chunk_size, 1);
        tmp_right_1d = zeros(total_chunk_size, 1);
        tmp_1d_percentile_mask = ones(total_chunk_size, 1);

        % multiply tmp_left and tmp_right by airvol to yield final molecule count (unit conversion necessary to sum to 99%)

        % read the data into left_1d and right_1d, multiplying by airvol ... to yield final molecule count
        current_pointer = 1;
        for I = 1:IM
        for J = 1:JM
            chunk_size = chunk_diff(I,J);
            tmp_left_1d(current_pointer:current_pointer+chunk_size-1)  = tmp_left(I,J,read_layer_ks(I,J):read_layer_ke(I,J)) .* airvol_3D(I,J,read_layer_ks(I,J):read_layer_ke(I,J));
            tmp_right_1d(current_pointer:current_pointer+chunk_size-1) = tmp_right(I,J,read_layer_ks(I,J):read_layer_ke(I,J)) .* airvol_3D(I,J,read_layer_ks(I,J):read_layer_ke(I,J));

            current_pointer = current_pointer + chunk_size;
        end
        end

        % here comes the tricky part.
        % Shen et al., 2022 states: the sum is over the Qi ordered grid boxes that account for 99% of the total mass of species i in the boundary layer (surface to 2km) ...
        % sorting is inevitable, so do the sorting
        [tmp_left_1d_sort, sort_idx] = sort(tmp_left_1d);
        tmp_right_1d_sort = tmp_right_1d(sort_idx); % need to match by sort_idx
        sum_target = 0.01 * sum(tmp_left_1d);
        sum_current = 0.0;

        % locate elements, marking as 0 to unmask while sum_target is not attained
        for I = 1:total_chunk_size
            sum_current = sum_current + tmp_left_1d_sort(I);
            tmp_1d_percentile_mask(I) = 0;
            if sum_current >= sum_target
                % fprintf("%s excluded bottom %d/%d grid boxes\n", rawSpcName, I, total_chunk_size);
                break
            end
        end
        tmp_1d_percentile_mask = logical(tmp_1d_percentile_mask);

        % also make a pass against > mask_threshold (10 molec cm-3)
        tmp_1d_percentile_mask = tmp_1d_percentile_mask & tmp_left_1d_sort > mask_threshold;

        % now we have the mask, calculate the RRMS ...
        % (median not supported here)
        rmse = sqrt(mean((tmp_right_1d_sort(tmp_1d_percentile_mask) - tmp_left_1d_sort(tmp_1d_percentile_mask)) .^ 2));
        rrms = sqrt(mean(((tmp_right_1d_sort(tmp_1d_percentile_mask) - tmp_left_1d_sort(tmp_1d_percentile_mask))./tmp_left_1d_sort(tmp_1d_percentile_mask)) .^ 2)) * 100;

        perspc_rrms(i) = rrms;
        perspc_rmse(i) = rmse;
        perspc_rrms_byName(rawSpcName) = rrms;

        if exist('external_log_rrms_perspc', 'var')
            external_log_rrms_perspc(rawSpcName)  = [external_log_rrms_perspc(rawSpcName), rrms];
            external_log_rmse_perspc(rawSpcName)  = [external_log_rmse_perspc(rawSpcName), rmse];
            external_log_val_left_mean_perspc(rawSpcName) = [external_log_val_left_mean_perspc(rawSpcName), mean(tmp_left_1d_sort(tmp_1d_percentile_mask))];
            external_log_val_right_mean_perspc(rawSpcName) = [external_log_val_right_mean_perspc(rawSpcName), mean(tmp_right_1d_sort(tmp_1d_percentile_mask))];
        end

        fprintf("=> %s RMSE = %.4e, RRMS = %.4f%s\n", rawSpcName, rmse, rrms, '%');
    % --------------------- SPLIT MODE ---------------------
    elseif strcmp(mode, "split")
        rmses = zeros(LM, 1);
        rrmss = zeros(LM, 1);
        for L = 1:LM
            tmp_left_1d = reshape(tmp_left(ips:ipe,jps:jpe,L), 1, IM*JM);
            tmp_right_1d = reshape(tmp_right(ips:ipe,jps:jpe,L), 1, IM*JM);

            mask = tmp_left_1d > mask_threshold;

            rmses(L) = sqrt(mean((tmp_right_1d(mask) - tmp_left_1d(mask)) .^ 2));
            rrmss(L) = sqrt(mean(((tmp_right_1d(mask) - tmp_left_1d(mask))./tmp_left_1d(mask)) .^ 2)) * 100;
        end

        % eventually check for do_flip_vert for plotting ...
        loglog(rrmss, levs * 1000, 'LineWidth', 2);
        xlabel("Relative Root Mean Squared Error (RRMS) [%]", 'FontSize', 14);
        ylabel("Pressure [hPa]", 'FontSize', 14);
        xlim([1e-5 100]);
        title(sprintf("Species = %s", rawSpcName), 'FontSize', 14);
        subtitle(sprintf("%s vs %s", filename_left, filename_right), 'Interpreter', 'none', 'FontSize', 10);
        set(gca, 'Ydir', 'reverse'); % plot pressure

        set(gcf, 'Renderer', 'painters', 'Position', [90 0 800 640]);

        % just save
        fname = sprintf("RRMS_byLevel_%s.png", rawSpcName);
        saveas(gcf, sprintf("out/%s/%s", casename, fname));
        close(gcf);
    elseif mode == "all"
        % for all, save into a 1-d array...
        % for first iteration, precompute size
        unified_chunk_size = IM*JM*LM;
        unified_chunk_size_2d = IM*JM;
        if i == 1
            % precompute size of arrays to prevent heap-allocation
            unified_3d_array_left = zeros(unified_chunk_size*(length(vars_in_3D)-1), 1);
            unified_3d_array_right = zeros(unified_chunk_size*(length(vars_in_3D)-1), 1);

            ptr = 1;
        else
            ptr = ptr + 1;
        end

        if mask_method == "relative"
            % load into chunks ...
            tmp_left_1d  = zeros(unified_chunk_size, 1);
            tmp_right_1d = zeros(unified_chunk_size, 1);

            for L = 1:LM
                tmp_left_2d = tmp_left(ips:ipe,jps:jpe,L);
                l_avg_thres = mask_threshold_rel_ratio * mean(tmp_left_2d(:));
                % l_avg_thres = 1e6;
                tmp_left_2d(tmp_left_2d <= l_avg_thres) = -999;

                tmp_right_2d = tmp_right(ips:ipe,jps:jpe,L);
                tmp_right_2d(tmp_left_2d <= l_avg_thres) = -999;

                % because a mask needs to be applied later, set -999
                tmp_left_1d(1+unified_chunk_size_2d*(L-1):unified_chunk_size_2d*L) = reshape(tmp_left_2d, 1, unified_chunk_size_2d);
                tmp_right_1d(1+unified_chunk_size_2d*(L-1):unified_chunk_size_2d*L) = reshape(tmp_right_2d, 1, unified_chunk_size_2d);

                % fprintf("%s - L = %d - Rel Thres %.6e\n", spiename, L, l_avg_thres);
            end

            unified_3d_array_left(1+unified_chunk_size*(ptr-1):unified_chunk_size*ptr) = tmp_left_1d(:);
            unified_3d_array_right(1+unified_chunk_size*(ptr-1):unified_chunk_size*ptr) = tmp_right_1d(:);
        else
            unified_3d_array_left(1+unified_chunk_size*(ptr-1):unified_chunk_size*ptr) = tmp_left(:);
            unified_3d_array_right(1+unified_chunk_size*(ptr-1):unified_chunk_size*ptr) = tmp_right(:);
        end

        if mod(i, 5) == 0
            fprintf("%d ", i)
        end
    elseif strcmp(mode, "map")
        % vertically integrated maps for RMSE (maybe split by thres_layer?)
        % make RMSE map (1:IM, 1:JM) for each and every vertical grid
        % box...
        if horiz_subset
            error("Mode == map cannot be used with horiz_subset");
        end

        % load into chunks ...
        tmp_left_1dv  = zeros(LM, 1);
        tmp_right_1dv = zeros(LM, 1);

        rrms_plt2d = zeros(IM, JM);

        % compute by-layer average for relative thresholding
        if strcmp(mask_method, "relative")
            fprintf("Using relative thresholding ... ");
            l_avg_thres = zeros(LM); % level average by species ...
            for L = 1:LM
                %fprintf("%d ", L);
                l_avg_thres(L) = mean(tmp_left(:,:,L), 'all') * mask_threshold_rel_ratio;

                % go through layers and remove those who are lower v. thres
                % since matlab does not have x(:,:,L)(x<y) construct
                % we need to do this the old way... very inefficient
                for I = 1:IM
                for J = 1:JM
                    if tmp_left(I,J,L) < l_avg_thres
                        tmp_left(I,J,L) = -999;
                        tmp_right(I,J,L) = -999;
                    end
                end
                end
            end
            fprintf("\n");
        end

        for I = 1:IM
        for J = 1:JM
            % still have to iterate through layers and vertically integ.
            % RRMS calculations
            % into rrms_plt2d
            tmp_left_1dv = squeeze(tmp_left(I,J,:));
            tmp_right_1dv = squeeze(tmp_right(I,J,:));
            mask = tmp_left_1dv > mask_threshold;

            if length(tmp_left_1dv(mask)) > 0
                rrms_plt2d(I,J) = sqrt(mean(((tmp_right_1dv(mask) - tmp_left_1dv(mask))./tmp_left_1dv(mask)) .^ 2)) * 100;
            else
                rrms_plt2d(I,J) = -1; % eliminated by mask
            end

            if isnan(rrms_plt2d(I,J))
                error("FATAL: NaN detected in loop ... check (1)")
            end
        end
        end

        if length(vars_in_3D) == 9
            % make subplots.
            % requires subplot_tight
            axg = subplot_tight(3,3,i);
        else
            % make one and save.
            axg = subplot(1,1,1);
        end

        hold on;
        m_proj('miller', ...
               'long',[-180.0 180.0], ...
               'lat', [-87.5 87.5]);

        m_pcolor(lons_2d, lats_2d, rrms_plt2d);
        cbr = colorbar;

        % check max val?
        if max(rrms_plt2d, [], 'all') < 1
            caxis([0, 1]);
        elseif max(rrms_plt2d, [], 'all') < 5
            caxis([0, 1]);
        elseif max(rrms_plt2d, [], 'all') < 10
            caxis([0, 10]);
        elseif max(rrms_plt2d, [], 'all') < 50
            caxis([0, 50]);
        elseif max(rrms_plt2d, [], 'all') < 100
            caxis([0, 100]);
        elseif max(rrms_plt2d, [], 'all') < 200
            caxis([0, 100]);
        elseif max(rrms_plt2d, [], 'all') < 500
            caxis([0, 500]);
        elseif max(rrms_plt2d, [], 'all') < 1000
            caxis([0, 1000]);
        elseif max(rrms_plt2d, [], 'all') < 2000
            caxis([0, 2000]);
        elseif max(rrms_plt2d, [], 'all') < 5000
            caxis([0, 5000]);
        else
            %caxis([0, 100]);
        end

        m_plotbndry(world_boundary_file, 'color', 'k', 'linewidth', 1.5);
        m_grid;
        colormap(axg, flip(othercolor('Spectral6', 10)));
        title(sprintf("RRMS [%s] - %s - %s", '%', rawSpcName, timestr), 'Interpreter', 'none');
        set(gcf, 'Renderer', 'painters', 'Position', [90 0 1500 1200]);

        if length(vars_in_3D) == 9
            fname = sprintf("RRMS_map9_%s_%s.png", casename, timestr);
        else
            fname = sprintf("RRMS_map_%s_%s_%s.png", casename, timestr, rawSpcName);
        end

        if length(vars_in_3D) ~= 9 || i == 9
            sgtitle(sprintf("RRMS Vertical [%s] (vs. %s) - %s", '%', casename, timestr), 'Interpreter', 'none');

            saveas(gcf, sprintf("out/%s", fname));
            close(gcf);
            fprintf("\nWritten to: out/%s\n", fname);
        end

    elseif strcmp(mode, "histogram")
        % separate surface - 100 hPa (layers 1:35) and above (geos-chem specific)
        thres_layer = 35;
        rrms_trop  = zeros(IM*JM*thres_layer, 1);
        rrms_strat = zeros(IM*JM*(LM-thres_layer), 1);
        mask_trop  = false(IM*JM*thres_layer, 1);
        mask_strat = false(IM*JM*(LM-thres_layer), 1);

        ptr = 1;


        for L = 1:LM
            tmp_left_1d = reshape(tmp_left(ips:ipe,jps:jpe,L), 1, IM*JM);
            tmp_right_1d = reshape(tmp_right(ips:ipe,jps:jpe,L), 1, IM*JM);

            if L == (thres_layer+1)
                ptr = 1;
            end

            if L <= thres_layer
                rrms_trop(ptr:ptr+IM*JM-1) = sqrt(((tmp_left_1d - tmp_right_1d) ./ tmp_left_1d) .^ 2) * 100;
                mask_trop(ptr:ptr+IM*JM-1) = tmp_left_1d > mask_threshold;
                if mask_method == "relative"
                    l_avg_thres = mask_threshold_rel_ratio * mean(tmp_left_1d(:));
                    mask_trop(ptr:ptr+IM*JM-1) = tmp_left_1d > mask_threshold & tmp_left_1d > l_avg_thres;
                end
                % fprintf("fill %d %d\n", ptr, ptr+IM*JM-1)
            else
                rrms_strat(ptr:ptr+IM*JM-1) = sqrt(((tmp_left_1d - tmp_right_1d) ./ tmp_left_1d) .^ 2) * 100;
                mask_strat(ptr:ptr+IM*JM-1) = tmp_left_1d > mask_threshold;
                if mask_method == "relative"
                    l_avg_thres = mask_threshold_rel_ratio * mean(tmp_left_1d(:));
                    mask_strat(ptr:ptr+IM*JM-1) = tmp_left_1d > mask_threshold & tmp_left_1d > l_avg_thres;
                end
            end

            ptr = ptr + IM*JM;
        end

        if length(vars_in_3D) == 9
            % make histogram plot now
            subplot(3,3,i);

            % calculate overall rrms with mask
            overall_rrms = sqrt((sum(rrms_trop(mask_trop) .^ 2) + sum(rrms_strat(mask_strat) .^ 2))/(IM*JM*LM));

            % manually lump (non-log)
            if histogram_scale == "regular"
                rrms_trop(rrms_trop > 200) = 200;
                rrms_strat(rrms_strat > 200) = 200;
            end

            % do a log10
            if histogram_scale == "log"
                [~,edges] = histcounts(log10(rrms_trop(:,1)), 100);
            end

            hold off;
            if histogram_scale == "regular"
                h1 = histogram(rrms_trop(:,1), 100, 'BinLimits', [0, 200]);
            end
            if histogram_scale == "log"
                h1 = histogram(rrms_trop(:,1), 10.^edges, 'EdgeAlpha', 0.2);
                xlim([10^-6 10^10]);
                % ylim([0 0.2*IM*JM*thres_layer]);
                [h1_maxCount, h1_whichBin] = max(h1.Values);
            end
            hold on;
            if histogram_scale == "regular"
                h2 = histogram(rrms_strat(:,1), 100, 'BinLimits', [0, 200]);
            end
            if histogram_scale == "log"
                h2 = histogram(rrms_strat(:,1), 10.^edges, 'EdgeAlpha', 0.2);
                % xlim([10^-6 10^10]);
                [h2_maxCount, h2_whichBin] = max(h2.Values);
            end
            legend("below 100 hPa", "above 100 hPa");
            title(sprintf("%s (%.2f%s)", rawSpcName, overall_rrms, '%'), 'Interpreter', 'none')

            if histogram_scale == "regular"
                subtitle("Values > 200 are lumped to 200")
            end
            if histogram_scale == "log"
                set(gca, 'xscale','log')
                subtitle(sprintf("trop: max ~ %2.2e c = %02d\nstrat: max ~ %2.2e c = %02d", 0.5*(h1.BinEdges(h1_whichBin)+h1.BinEdges(h1_whichBin+1)), h1_maxCount, ...
                                                                                            0.5*(h2.BinEdges(h2_whichBin)+h2.BinEdges(h2_whichBin+1)), h2_maxCount), 'FontName', 'DejaVu Sans Mono');

                if exist('external_log_histogram_peak_trop', 'var')
                    external_log_histogram_peak_trop(rawSpcName)  = [external_log_histogram_peak_trop(rawSpcName),  0.5*(h1.BinEdges(h1_whichBin)+h1.BinEdges(h1_whichBin+1))];
                    external_log_histogram_peak_strat(rawSpcName) = [external_log_histogram_peak_strat(rawSpcName), 0.5*(h2.BinEdges(h2_whichBin)+h2.BinEdges(h2_whichBin+1))];
                end
            end

            if i == 9
                set(gcf, 'Renderer', 'painters', 'Position', [90 0 1500 1200]);

                % just save
                if horiz_subset
                    extra = sprintf("I_%d-%d_J_%d-%d", ips, ipe, jps, jpe);
                    extra_nice = sprintf("I:%d-%d J:%d-%d", ips, ipe, jps, jpe);
                else
                    extra = "";
                    extra_nice = "Whole grid";
                end
                sgtitle(sprintf("RRMS Histogram [%s] %s (vs. %s) - %s", '%', extra_nice, casename, timestr), 'Interpreter', 'none');

                fname = sprintf("RRMS_histogram9_%s_%s%s.png", casename, timestr, extra);
                saveas(gcf, sprintf("out/%s", fname));
                close(gcf);

                fprintf("\nWritten to: out/%s\n", fname);
            end
        else
            error("TODO: mode == histogram has to specify exactly 9 species for now");
        end
    end
end

if mode == "all"
    % compute the total RMSE, RRMS
    %
    % according to Santillana et al., 2010; Shen et al., 2020, the mask has
    % a threshold of 1e6
    if mask_method == "relative"
        mask = unified_3d_array_left > mask_threshold; % or 0
    else
        mask = unified_3d_array_left > mask_threshold;
    end

    if use_median_not_mean == 1
        rmse = sqrt(median((unified_3d_array_right(mask) - unified_3d_array_left(mask)) .^ 2));
        rrms = sqrt(median(((unified_3d_array_right(mask) - unified_3d_array_left(mask))./unified_3d_array_left(mask)) .^ 2)) * 100;
    else
        rmse = sqrt(mean((unified_3d_array_right(mask) - unified_3d_array_left(mask)) .^ 2));
        rrms = sqrt(mean(((unified_3d_array_right(mask) - unified_3d_array_left(mask))./unified_3d_array_left(mask)) .^ 2)) * 100;
    end

    fprintf("\n\n=> [All Species, All Levels] RMSE = %.8e, RRMS = %.4f%s\n", rmse, rrms, '%');

    if exist('external_log_rrms_all', 'var')
        external_log_rrms_all = [external_log_rrms_all, rrms];
    end
end

if (strcmp(mode, "perspc") || strcmp(mode, "perspc2")) && ~exist('external_log_rmse_perspc', 'var')
    % sort ...
    [perspc_rrms_sort, sort_idx] = sort(perspc_rrms, 'descend');
    vars_in_3D_sortrrms = vars_in_3D(sort_idx);
    perspc_rmse_sortrrms = perspc_rmse(sort_idx);

    fprintf("======= Sorted List of RRMS (descending) =======\n");
    for i = 1:length(vars_in_3D_sortrrms)
        fprintf("=> %s RMSE = %.8e, RRMS = %.4f%s\n", vars_in_3D_sortrrms(i), perspc_rmse_sortrrms(i), perspc_rrms_sort(i), '%');
    end


    fprintf("Comparing:\nLeft_Model: %s\nRight_Model: %s\n", filename_left, filename_right);
    if strcmp(mode, "perspc")
        fprintf("perspc_flag: %d (5=500hPa,6=sfc,7=strat,8=trop)\n", perspc_flag);
    elseif strcmp(mode, "perspc2")
        fprintf("perspc2_subset: %s\n", perspc2_subset);
    end
    fprintf("======= Benchmark category results =======\n");
    aerosols_dust = nanmean([perspc_rrms_byName("DST1"), perspc_rrms_byName("DST2"), perspc_rrms_byName("DST3"), perspc_rrms_byName("DST4")]);
    aerosols_inorganic = nanmean([perspc_rrms_byName("NH4"), perspc_rrms_byName("NIT"), perspc_rrms_byName("SO4")]);
    aerosols_ocbc = nanmean([perspc_rrms_byName("BCPI"), perspc_rrms_byName("BCPO"), perspc_rrms_byName("OCPI"), perspc_rrms_byName("OCPO")]);
    aerosols_seasalt = nanmean([perspc_rrms_byName("AERI"), perspc_rrms_byName("BrSALA"), perspc_rrms_byName("BrSALC"), perspc_rrms_byName("ISALA"), perspc_rrms_byName("ISALC"), perspc_rrms_byName("NITs"), perspc_rrms_byName("SALA"), perspc_rrms_byName("SALAAL"), perspc_rrms_byName("SALACL"), perspc_rrms_byName("SALC"), perspc_rrms_byName("SALCAL"), perspc_rrms_byName("SALCCL"), perspc_rrms_byName("SO4s")]);
    aerosols = nanmean([perspc_rrms_byName("DST1"), perspc_rrms_byName("DST2"), perspc_rrms_byName("DST3"), perspc_rrms_byName("DST4"), perspc_rrms_byName("NH4"), perspc_rrms_byName("NIT"), perspc_rrms_byName("SO4"), perspc_rrms_byName("BCPI"), perspc_rrms_byName("BCPO"), perspc_rrms_byName("OCPI"), perspc_rrms_byName("OCPO"), perspc_rrms_byName("AERI"), perspc_rrms_byName("BrSALA"), perspc_rrms_byName("BrSALC"), perspc_rrms_byName("ISALA"), perspc_rrms_byName("ISALC"), perspc_rrms_byName("NITs"), perspc_rrms_byName("SALA"), perspc_rrms_byName("SALAAL"), perspc_rrms_byName("SALACL"), perspc_rrms_byName("SALC"), perspc_rrms_byName("SALCAL"), perspc_rrms_byName("SALCCL"), perspc_rrms_byName("SO4s"), perspc_rrms_byName("pFe")]);
    aerosols_metals = nanmean([perspc_rrms_byName("pFe")]);
    bromine = nanmean([perspc_rrms_byName("Br"), perspc_rrms_byName("Br2"), perspc_rrms_byName("BrCl"), perspc_rrms_byName("BrNO2"), perspc_rrms_byName("BrNO3"), perspc_rrms_byName("BrO"), perspc_rrms_byName("CH3Br"), perspc_rrms_byName("CH2Br2"), perspc_rrms_byName("CHBr3"), perspc_rrms_byName("HOBr"), perspc_rrms_byName("HBr")]);
    chlorine = nanmean([perspc_rrms_byName("Cl"), perspc_rrms_byName("ClO"), perspc_rrms_byName("Cl2"), perspc_rrms_byName("Cl2O2"), perspc_rrms_byName("ClOO"), perspc_rrms_byName("ClNO2"), perspc_rrms_byName("ClNO3"), perspc_rrms_byName("CCl4"), perspc_rrms_byName("CH3Cl"), perspc_rrms_byName("CH2Cl2"), perspc_rrms_byName("CH3CCl3"), perspc_rrms_byName("CHCl3"), perspc_rrms_byName("HOCl"), perspc_rrms_byName("HCl"), perspc_rrms_byName("OClO")]);
    iodine = nanmean([perspc_rrms_byName("I"), perspc_rrms_byName("I2"), perspc_rrms_byName("IBr"), perspc_rrms_byName("ICl"), perspc_rrms_byName("IO"), perspc_rrms_byName("IONO"), perspc_rrms_byName("IONO2"), perspc_rrms_byName("CH3I"), perspc_rrms_byName("CH2I2"), perspc_rrms_byName("CH2ICl"), perspc_rrms_byName("CH2IBr"), perspc_rrms_byName("HI"), perspc_rrms_byName("HOI"), perspc_rrms_byName("OIO")]);
    nitrogen = nanmean([perspc_rrms_byName("HNO2"), perspc_rrms_byName("HNO4"), perspc_rrms_byName("MPAN"), perspc_rrms_byName("HNO3T"), perspc_rrms_byName("NO"), perspc_rrms_byName("NO2"), perspc_rrms_byName("NO3"), perspc_rrms_byName("N2O5"), perspc_rrms_byName("MPN"), perspc_rrms_byName("PAN"), perspc_rrms_byName("PPN"), perspc_rrms_byName("N2O"), perspc_rrms_byName("NH3T"), perspc_rrms_byName("MENO3"), perspc_rrms_byName("ETNO3"), perspc_rrms_byName("IPRNO3"), perspc_rrms_byName("NPRNO3")]);
    oxidants = nanmean([perspc_rrms_byName("O3"), perspc_rrms_byName("CO"), perspc_rrms_byName("OH"), perspc_rrms_byName("NO2"), perspc_rrms_byName("NO"), perspc_rrms_byName("NO3")]); % HNO2, HNO4, INO, IONO, IONO2, N, OLND, OLNN
    primary_org_alcohols = nanmean([perspc_rrms_byName("EOH"), perspc_rrms_byName("MOH")]);
    primary_org_biogenics = nanmean([perspc_rrms_byName("ISOP"), perspc_rrms_byName("MTPA"), perspc_rrms_byName("MTPO"), perspc_rrms_byName("LIMO")]);
    primary_org_hcs = nanmean([perspc_rrms_byName("ALK4"), perspc_rrms_byName("BENZ"), perspc_rrms_byName("CH4"), perspc_rrms_byName("C2H6"), perspc_rrms_byName("C3H8"), perspc_rrms_byName("PRPE"), perspc_rrms_byName("TOLU"), perspc_rrms_byName("XYLE")]);
    primary_org = nanmean([perspc_rrms_byName("EOH"), perspc_rrms_byName("MOH"), perspc_rrms_byName("ISOP"), perspc_rrms_byName("MTPA"), perspc_rrms_byName("MTPO"), perspc_rrms_byName("LIMO"), perspc_rrms_byName("ALK4"), perspc_rrms_byName("BENZ"), perspc_rrms_byName("CH4"), perspc_rrms_byName("C2H6"), perspc_rrms_byName("C3H8"), perspc_rrms_byName("PRPE"), perspc_rrms_byName("TOLU"), perspc_rrms_byName("XYLE")]);
    roy = nanmean([perspc_rrms_byName("H2O2"), perspc_rrms_byName("H"), perspc_rrms_byName("H2"), perspc_rrms_byName("H2O"), perspc_rrms_byName("HO2"), perspc_rrms_byName("O1D"), perspc_rrms_byName("OH")]); % , perspc_rrms_byName("RO2")); RO2 is lumped
    secondary_organics = nanmean([perspc_rrms_byName("ACTA"), perspc_rrms_byName("ALD2"), perspc_rrms_byName("CH2O"), perspc_rrms_byName("MACR"), perspc_rrms_byName("HPALD1"), perspc_rrms_byName("HPALD2"), perspc_rrms_byName("HPALD3"), perspc_rrms_byName("HPALD4"), perspc_rrms_byName("IEPOXA"), perspc_rrms_byName("IEPOXB"), perspc_rrms_byName("IEPOXD"), perspc_rrms_byName("ACET"), perspc_rrms_byName("HAC"), perspc_rrms_byName("MEK"), perspc_rrms_byName("MVK"), perspc_rrms_byName("IHN1"), perspc_rrms_byName("IHN2"), perspc_rrms_byName("IHN3"), perspc_rrms_byName("IHN4"), perspc_rrms_byName("GLYC"), perspc_rrms_byName("GLYX"), perspc_rrms_byName("HCOOH"), perspc_rrms_byName("MAP"), perspc_rrms_byName("RCHO"), perspc_rrms_byName("MP")]);
    soa = nanmean([perspc_rrms_byName("SOAP"), perspc_rrms_byName("SOAS")]);
    sulfur = nanmean([perspc_rrms_byName("DMS"), perspc_rrms_byName("MSA"), perspc_rrms_byName("OCS"), perspc_rrms_byName("SO2"), perspc_rrms_byName("SO4")]);

    fprintf("=> Oxidants: %.2f%s\n", oxidants, '%');
    fprintf("=> ROy: %.2f%s\n", roy, '%');
    fprintf("=> Aerosols: %.2f%s\n", aerosols, '%');
    fprintf("   - Dust: %.2f%s\n", aerosols_dust, '%');
    fprintf("   - Inorganic: %.2f%s\n", aerosols_inorganic, '%');
    fprintf("   - OC/BC: %.2f%s\n", aerosols_ocbc, '%');
    fprintf("   - SeaSalt: %.2f%s\n", aerosols_seasalt, '%');
    fprintf("   - Metals: %.2f%s\n", aerosols_metals, '%');
    fprintf("=> Nitrogen: %.2f%s\n", nitrogen, '%');
    fprintf("=> Primary Organics: %.2f%s\n", primary_org, '%');
    fprintf("   - Alcohols: %.2f%s\n", primary_org_alcohols, '%');
    fprintf("   - Biogenics: %.2f%s\n", primary_org_biogenics, '%');
    fprintf("   - HCs: %.2f%s\n", primary_org_hcs, '%');
    fprintf("=> Secondary Organics: %.2f%s\n", secondary_organics, '%');
    fprintf("=> SOA: %.2f%s\n", soa, '%');
    fprintf("=> Sulfur: %.2f%s\n", sulfur, '%');
    fprintf("=> Br: %.2f%s\n", bromine, '%');
    fprintf("=> Cl: %.2f%s\n", chlorine, '%');
    fprintf("=> I: %.2f%s\n", iodine, '%');

end

if isautoallflag == true
    % calculate the mean and median of all spec rrms.
    %[perspc_rrms_sort, sort_idx] = sort(perspc_rrms, 'descend');
    %vars_in_3D_sortrrms = vars_in_3D(sort_idx);

    perspc_rrms(isnan(perspc_rrms)) = []; % remove NaN

    fprintf("\n==> All spc RRMS mean: %.4f, 5th: %.2f, 25th: %.2f, median: %.2f, 75th: %.2f, 95th: %.2f\n", mean(perspc_rrms), prctile(perspc_rrms,5), prctile(perspc_rrms,25), median(perspc_rrms), prctile(perspc_rrms,75), prctile(perspc_rrms,95));

    if exist('external_log_rrms_allspc_mean', 'var')
        external_log_rrms_allspc_mean = [external_log_rrms_allspc_mean, mean(perspc_rrms)];
        external_log_rrms_allspc_median = [external_log_rrms_allspc_median, median(perspc_rrms)];

        external_log_rrms_allspc_prc05 = [external_log_rrms_allspc_prc05, prctile(perspc_rrms, 5)];
        external_log_rrms_allspc_prc25 = [external_log_rrms_allspc_prc25, prctile(perspc_rrms, 25)];
        external_log_rrms_allspc_prc75 = [external_log_rrms_allspc_prc75, prctile(perspc_rrms, 75)];
        external_log_rrms_allspc_prc95 = [external_log_rrms_allspc_prc95, prctile(perspc_rrms, 95)];
    end
end

% fprintf("\nEnd Script\n");


% Changelog:
% 2022.04.20 - Add perspc2, using Shen et al. 2022 filtering. Very slow, but gets species that sum up to 99% of the mass
% 2022.04.15 - Add benchmark species
% 2022.04.03 - Map for one level only
% 2022.02.25 - Merge NH3T = NH3+NH4, HNO3T = HNO3+NIT; for perspc-autoall, show all spc rrms mean and median
% 2022.02.16 - Support Median for "all" type -- use median of entire 1-D collapsed array?
% 2022.01.31 - Strat|Trop (L = 35) split for calculating perspc rrms/rmse.
% 2022.01.13 - Save mean conc corresponding to calculated rmse as well for test.
% 2022.01.12 - Save rmse as well. Change diff to rrms for clarity.
% 2022.01.11 - Add common spc; allow saving all/perspc into forced external ts.
% 2021.12.09 - Now sort RRMS for perspc.
% 2021.12.07 - Find histogram peak
% 2021.12.06 - Histogram saving; multiple time slice in file stub (no
% support); external forcing of date stamp
% 2021.11.21 - Now
% 2021.11.16 - Initial version
% Laura Palmer - Bastille
% Cambridge MA

%--- tip
% to force timestr externally:
% clear external_log_histogram_peak_trop;
% clear external_log_rrms_all;
% clear external_log_rrms_perspc;
% timestrs1 = arrayfun(@(t){sprintf("20190101_%02d00", t)}, 3:3:21);
% timestrs = [timestrs1];
% for i = 1:12
% for j = 1:31
%     if i == 1 && j == 1
%         continue
%     end
%     if i == 2 && j > 28
%         continue
%     end
%     if i == 4 && j > 30
%         continue
%     end
%     if i == 6 && j > 30
%         continue
%     end
%     if i == 9 && j > 30
%         continue
%     end
%     if i == 11 && j > 30
%         continue
%     end
%     timestrs_tmp = arrayfun(@(t){sprintf("2019%02d%02d_%02d00", i, j, t)}, 0:3:21);
%     timestrs = [timestrs, timestrs_tmp];
% end
% end
% timestrs_final = [timestrs, "20200101_0000"];
% for i = 1:numel(timestrs_final)
%     external_timestr = timestrs_final{i}; run("Metrics_Diff_RMSE.m")
% end

% save("outs_longrun_1yr_all.mat")
%--- /end

