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
%
% Reads raw netCDF files for convenience
%
% (c) 2019-2021 Haipeng Lin <jimmie.lin@gmail.com>
%
% Version: 2021.12.07
% Started: 2021.11.16
% See changelog at end of script
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% specify data on left and right
%
% netCDF file names (standard model output)
%
% RMSE percentages (RRMS) are calculated based on filename_left.
casename = 'new_diag_t100_prs_1day_multi';
base_casename = 'new_diag_baseline_1day_multi';

casename = 'v3_t100_prs';
base_casename = 'v3_baseline';

if exist('external_timestr', 'var')
    timestr = external_timestr;
else
    timestr  = '20190708_0000';
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
%vars_in_3D = ["SpeciesConcMND_Br", "SpeciesConcMND_BrO", "SpeciesConcMND_HBr", ...
%              "SpeciesConcMND_Cl", "SpeciesConcMND_ClO", "SpeciesConcMND_HCl", ...
%              "SpeciesConcMND_O3", "SpeciesConcMND_NO2", "SpeciesConcMND_CO"];

vars_in_3D = ["SpeciesConcMND_Br", "SpeciesConcMND_ALD2"]; % testing only

% mode: unified|split
%
% all:     compute RMSE over the entire vertical for ALL 3-D SpeciesConcMND_
% all-lev: compute RMSE over each level individually for ALL 3-D SpeciesConcMND_
% perspc:  compute RMSE over the entire vertical for some vars
%  perspc-autoall: all, but no need to fill vars_in_3D yourself
% split:   compute RMSE over each level individually plotting end result, for some vars
% histogram: RRMS over every grid box and plot histogram
% map:     make a m_map (boundary package required) for RRMS, vertically integrated
mode = "perspc-autoall";

% histogram_scale: log | regular (if mode == 'histogram')
histogram_scale = "log";

% mask_method: relative|absolute
%
% if absolute, any value_left < mask_threshold will be ignored
% if relative, any value_left < rel_ratio * vertical layer mean will be ignored
%   rel_ratio usually chosen 0.1 (ignore lower than 10% * avg) or 0.01
%            , AND any value_left < mask_threshold will be ignored
mask_method = "relative";

% mask_threshold (if mask_method == 'absolute')
% set to 1e6 molec/cm3 per Santillana et al., 2010; Shen et al., 2020
%
% since data is in v/v dry, convert ppv to molec/cm3 roughly using
% 1 atm ~ 2.46e19 molec/cm3
%
%
% used to avoid division by zero
mask_threshold = 1e4;
mask_threshold_rel_ratio = 0.1;

% IM, JM subsetting
% do any subsetting of data? specify range if so.
%
% ignored for world plots obviously
%
% warning: no effort is made to check IM_, JM_range validity
horiz_subset = false;
IM_range = 13:17;
JM_range = 23:29;

% world_boundary_file: '/usr/local/MATLAB/R2021b/toolbox/boundary/world'
% get this from the boundary toolbox.
world_boundary_file = '/usr/local/MATLAB/R2021b/toolbox/boundary/world';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% NO USER CONFIGURABLE CODE BELOW %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf("Comparing:\nLeft_Model: %s\nRight_Model: %s\n", filename_left, filename_right);

% if all species, replace vars_in_3D with all the SpeciesConcMND_
%
if strcmp(mode, "all") || strcmp(mode, "perspc-autoall")
    vars_in_3D = ["dummy"]; % because otherwise we cannot define an all-string list...
    vars_raw_nc = {ncinfo(filename_left).Variables.Name};
    for i = 1:length(vars_raw_nc)
        if contains(char(vars_raw_nc(i)), "SpeciesConcMND_")
            % fprintf("check: %s\n", char(vars_raw_nc(i)))
            %if contains(char(vars_raw_nc(i)), "SpeciesConcMND_HNO3") || contains(char(vars_raw_nc(i)), "SpeciesConcMND_NH3")
            %    continue
            %end

            vars_in_3D(end+1) = sprintf("%s", char(vars_raw_nc(i)));
        end
    end

    if strcmp(mode, "perspc-autoall")
        mode = "perspc";
    end
end

% sort the species list in order
if ~strcmp(mode, "histogram")
    vars_in_3D = sort(vars_in_3D);
end

% if perspc, store rrms per species so we can sort later
if strcmp(mode, "perspc")
    perspc_rrms = zeros(length(vars_in_3D), 1);
    perspc_rmse = zeros(length(vars_in_3D), 1);
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
    tmp_left  = ncread(filename_left,  spiename);
    tmp_right = ncread(filename_right, spiename);
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
            tmp_left_1d  = zeros(unified_chunk_size*LM, 1);
            tmp_right_1d = zeros(unified_chunk_size*LM, 1);
    
            for L = 1:LM
                tmp_left_2d = tmp_left(ips:ipe,jps:jpe,L);
                l_avg_thres = mask_threshold_rel_ratio * mean(tmp_left_2d(:));
                % l_avg_thres = 1e6;
                tmp_left_2d(tmp_left_2d <= l_avg_thres) = -999;

                tmp_right_2d = tmp_right(ips:ipe,jps:jpe,L);
                tmp_right_2d(tmp_left_2d <= l_avg_thres) = -999;

                % because a mask needs to be applied later, set -999
                tmp_left_1d(1+unified_chunk_size*(L-1):unified_chunk_size*L) = reshape(tmp_left_2d, 1, unified_chunk_size);
                tmp_right_1d(1+unified_chunk_size*(L-1):unified_chunk_size*L) = reshape(tmp_right_2d, 1, unified_chunk_size);

                % fprintf("%s - L = %d - Rel Thres %.6e\n", spiename, L, l_avg_thres);
            end
        
            mask = tmp_left_1d > mask_threshold;

            rmse = sqrt(mean((tmp_right_1d(mask) - tmp_left_1d(mask)) .^ 2));
            rrms = sqrt(mean(((tmp_right_1d(mask) - tmp_left_1d(mask))./tmp_left_1d(mask)) .^ 2)) * 100;

            if isempty(tmp_left_1d(mask))
                fprintf("Warning: %s has no valid elements after threshold; before threshold min = %e\n", spiename, min(tmp_left_1d(:)));
                rmse = 0; rrms = 0;
            end
        else
            tmp_left_1d = tmp_left(:);
            tmp_right_1d = tmp_right(:);
        
            mask = tmp_left_1d > mask_threshold;

            rmse = sqrt(mean((tmp_right_1d(mask) - tmp_left_1d(mask)) .^ 2));
            rrms = sqrt(mean(((tmp_right_1d(mask) - tmp_left_1d(mask))./tmp_left_1d(mask)) .^ 2)) * 100;
        end

        perspc_rrms(i) = rrms;
        perspc_rmse(i) = rmse;

        fprintf("=> [Per-Species] %s RMSE = %.8e, RRMS = %.4f%s (minval=%e)\n", spiename, rmse, rrms, '%', min(tmp_left(:)));
    elseif mode == "split"
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
        if max(rrms_plt2d, [], 'all') < 50
            caxis([0, 50]);
        elseif max(rrms_plt2d, [], 'all') < 100
            caxis([0, 100]);
        elseif max(rrms_plt2d, [], 'all') < 200
            caxis([0, 200]);
        else
            caxis([0, 500]);
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

    rmse = sqrt(mean((unified_3d_array_right(mask) - unified_3d_array_left(mask)) .^ 2));
    rrms = sqrt(mean(((unified_3d_array_right(mask) - unified_3d_array_left(mask))./unified_3d_array_left(mask)) .^ 2)) * 100;

    fprintf("\n\n=> [All Species, All Levels] RMSE = %.8e, RRMS = %.4f%s\n", rmse, rrms, '%');
end

if strcmp(mode, "perspc")
    % sort ...
    [perspc_rrms_sort, sort_idx] = sort(perspc_rrms, 'descend');
    vars_in_3D_sortrrms = vars_in_3D(sort_idx);
    perspc_rmse_sortrrms = perspc_rmse(sort_idx);

    fprintf("======= Sorted List of RRMS (descending) =======\n");
    for i = 1:length(vars_in_3D_sortrrms)
        fprintf("=> %s RMSE = %.8e, RRMS = %.4f%s\n", vars_in_3D_sortrrms(i), perspc_rmse_sortrrms(i), perspc_rrms_sort(i), '%');
    end
end

fprintf("\nEnd Script\n");


% Changelog:
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
%timestrs1 = arrayfun(@(t){sprintf("20190701_%02d00", t)}, 3:3:21);
%timestrs2 = arrayfun(@(t){sprintf("20190702_%02d00", t)}, 0:3:21);
%timestrs = [timestrs1, timestrs2, "20190703_0000"]

%for i = 1:numel(timestrs_final)
%    external_timestr = timestrs_final{i}; run("Metrics_Diff_RMSE.m")
%end
%--- /end
