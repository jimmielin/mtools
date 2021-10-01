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
% (c) 2019-2021 Haipeng Lin <jimmie.lin@gmail.com>
%
% Version: dev
% Started: 2021.09.26
% See changelog at end of script
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% specify data coordinates. provide the "load" handle to the rest of the
% routines
coords_CESM = load('coords_CESM.mat');
coords_left = coords_CESM;
coords_right= coords_CESM;

% test species first; later generalize to all species
fname_left  = "climo_CAMchemYearlyMonthlyAvg_2016.mat";
fname_right = "climo_CESM2GCYearlyMonthlyAvg_2016.mat";

left_name      = "CESM2-GC";
right_name     = "CAM-chem";

% specify data indexing properties
%
% fourth_idx: usually date slice
fourth_idx = 7; % date index

% species to compare.
spec_list = ["CH2O", "CO", "HO2", "ISOP", "NO", "NO2", "O3", "OH", "SO2"];
%spec_list = ["O3"];

%%%%%%%%%%%%%%%%%%%%%% NO USER CONFIGURABLE CODE BELOW %%%%%%%%%%%%%%%%%%%%

%% load data
file_left   = load(fname_left);
file_right  = load(fname_right);

%% process data (in separate section to save time for debug)
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

psfc_left = file_left.PSFC(:,:,fourth_idx);
psfc_right = file_right.PSFC(:,:,fourth_idx);
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

% compute PEDGE
pedge_left = pedge_calc(coords_left.hyai, coords_left.hybi, psfc_left);
pedge_right = pedge_calc(coords_right.hyai, coords_right.hybi, psfc_right);

% compute DELP_DRY, which is just dry pedge deltas
delp_dry_left = delp_dry_calc(pedge_left);
delp_dry_right = delp_dry_calc(pedge_right);

% compute grid box mass for each ... using area_m2, delp_dry
AD_left = ad_calc(delp_dry_left, area_m2_left);
AD_right = ad_calc(delp_dry_right, area_m2_right);

% species specific code below. do loop
% sum species by multiplying [kg] in grid box by MMR
LM_left  = size(AD_left, 3);
LM_right = size(AD_right, 3);

% convert MMR to VMR
% VMR = 28.9644 / molar_mass * MMR
% MMR = molar_mass / 28.9644 * VMR
% a hack for a molar mass table ...
molarTable = containers.Map;
molarTable('CH2O') = 30.031;
molarTable('CO') = 28.01;
molarTable('HO2') = 33.01;
molarTable('ISOP') = 68.12;
molarTable('NO') = 30.01;
molarTable('NO2') = 46.0055;
molarTable('O3') = 48.00;
molarTable('OH') = 17.008;
molarTable('SO2') = 64.066;

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
    
    VMRtoMMR = molarTable(spec_list(idx)) / 28.9644;
    
    for L = 1:LM_left
        % sum by levels ..                 vv note this is element-wise mult.
        by_layer = sum(sum(AD_left(:,:,L) .* spec_current_left(:,:,L,fourth_idx) * VMRtoMMR));
        Msum_l = Msum_l + by_layer;
        %fprintf("=> [L] spc %s, L = %d, AD = %4.4e, v/v = %8.6e, m = %8.6f Gg\n", spec_list(idx), L, AD_left(25,25,L), spec_current_left(25,25,L,fourth_idx), by_layer/1e6);
    end

    for L = 1:LM_right
        % sum by levels ..                  vv note this is element-wise mult.
        by_layer = sum(sum(AD_right(:,:,L) .* spec_current_right(:,:,L,fourth_idx) * VMRtoMMR));
        Msum_r = Msum_r + by_layer;
        %fprintf("=> [R] spc %s, L = %d, v/v = %8.6e, m = %8.6f Gg\n", spec_list(idx), L, vv, by_layer/1e6);
    end

    % fprintf("Column sums L = %.6f Gg, R = %.6f Gg\n", Msum_l/1e6, Msum_r/1e6)
    fprintf("%-10s             %-8.6f \t\t%-8.6f \t\t%+02.3f\n", spec_list(idx), Msum_l/1e6, Msum_r/1e6, (Msum_l-Msum_r)/Msum_r)
end

% calculate ozone column value and plot (against climatology?)
% this is an added bonus IF the O3 field is present
if isfield(file_left, 'O3') && isfield(file_right, 'O3')
    fprintf("=> Also showing ozone column plot...\n");
    
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
        O3_col_l = O3_col_l + AD_left(:,:,L) .* spec_current_left(:,:,L,fourth_idx) ./ area_m2_left(:,:) * 46698 * 48 / 28.9644;
    end
    
    for L = 1:LM_right
        O3_col_r = O3_col_r + AD_right(:,:,L) .* spec_current_right(:,:,L,fourth_idx) ./ area_m2_right(:,:) * 46698 * 48 / 28.9644;
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
    axg = subplot_tight(1,2,1);

    hold on;
    m_proj('miller', ...
           'long',[-180.0 180.0], ...
           'lat', [-87.5 87.5]);

    m_pcolor(lons_left_sf, lats_left, O3_col_l);
    colorbar;
    caxis([150, 450]);

    m_plotbndry('/usr/local/MATLAB/R2021a/toolbox/boundary/world', ...
        'color', 'k', 'linewidth', 1.5)
    m_grid;
    colormap(axg, flip(othercolor('Spectral6', 12)));

    t = title(sprintf(left_name));

    % RIGHT
    axg = subplot_tight(1,2,2);

    hold on;
    m_proj('miller', ...
           'long',[-180.0 180.0], ...
           'lat', [-87.5 87.5]);

    m_pcolor(lons_right_sf, lats_right, O3_col_r);
    colorbar;
    caxis([150, 450]);

    m_plotbndry('/usr/local/MATLAB/R2021a/toolbox/boundary/world', ...
        'color', 'k', 'linewidth', 1.5)
    m_grid;
    colormap(axg, flip(othercolor('Spectral6', 12)));

    t = title(sprintf(right_name));
    sgtitle("O3 column [DU]");

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
