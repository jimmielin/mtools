%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%            _____            ______        %
% _______ _____  /_______________  /_______ %
% __  __ `__ \  __/  __ \  __ \_  /__  ___/ %
% _  / / / / / /_ / /_/ / /_/ /  / _(__  )  %
% /_/ /_/ /_/\__/ \____/\____//_/  /____/   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       "mtools" Research Toolkit           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Util_OMIMLSClimoASCII.m
%
% Reads OMI/MLS ozone column climatology data at 5x5 deg resolution
% (from https://acd-ext.gsfc.nasa.gov/Data_services/cloud_slice/)
% ASCII format files into matrices.
% (c) 2019-2021 Haipeng Lin <jimmie.lin@gmail.com>
%
% Version: 2021.10.01
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fid = fopen('sco_5by5.txt');

% create the matrix based on 5x5 resolution
delta = 5;
lons = -177.5:delta:180;
lats = -87.5:delta:90;
out_mat = zeros(size(lons, 2), size(lats, 2), 12); % 12-month climo.

% read line by line
tline = fgetl(fid);
curr_lat = -999;
while ischar(tline)
    % are we in a new latitude?
    if regexp(tline, "Latitude")
        % check if pos or neg
        l = tline;
        [tokens, matches] = regexp(l, "Latitude: (\d+)([NS])?-(\d+)([NS])?", 'tokens', 'match');
        if tokens{1}{4} == 'S'
            latcenter = (-1) * str2double(tokens{1}{3}) + delta/2;
        elseif tokens{1}{4} == 'N'
            latcenter = str2double(tokens{1}{3}) - delta/2;
        else
            error("Invalid tokens.")
        end
            
        fprintf("\n=> Reading Latitude line. Latitude starts: %f. Line is: %s\n", latcenter, l); 
        curr_lat = latcenter; % update current latitude
    elseif regexp(tline, "[WE0]-") % now try to match data up ...
        % read line split ... will be cell array
        linesplit = regexp(tline, '\s+', 'split');
        % iterate through each... each entry is one in 12 months
        for i = 1:size(linesplit, 2)
            currsplit = string(linesplit(i)); % R2016b+
            if regexp(currsplit, "[WE0]-")
                % match coords again
                [tokens, matches] = regexp(currsplit, "(\d+)([WE])?-(\d+)([WE])?", 'tokens', 'match');
                % 175W-180, 170W-175W, ... 0-5W, 0-5E, 5E-10E, ...
                % 175E-180E
                if size(tokens{1}{4}, 1) == 0 || tokens{1}{4} == 'W'
                    loncenter = (-1) * str2double(tokens{1}{3}) + delta/2;
                    fprintf("%fW ", loncenter);
                elseif tokens{1}{4} == 'E'
                    loncenter = str2double(tokens{1}{3}) - delta/2;
                    fprintf("%fE ", loncenter);
                else
                    error("Invalid tokens.")
                end
                
                % read the following 12 entries ...
                % linesplit(1, i+1:i+12)
                bymonthentry = cellfun(@str2num, linesplit(1, i+1:i+12));
                
                % locate the appropriate lons, lats index
                i_idx = find(lons==loncenter);
                j_idx = find(lats==curr_lat);
                out_mat(i_idx, j_idx, 1:12) = bymonthentry(1:12);
                fprintf("i_idx,j_idx=%d,%d ", i_idx, j_idx);
            end
        end
    end
    
    tline = fgetl(fid); % advance to next line
end

fclose(fid);
        
% pause
fprintf("\n");

%% plot test:
if size(lons, 1) == 1
    lons_new = zeros(size(lons, 2), size(lats, 2));
    lats_new = zeros(size(lons, 2), size(lats, 2));

    for j = 1:size(lats, 2)
        lons_new(:,j) = lons(:);
    end

    for i = 1:size(lons, 2)
        lats_new(i,:) = lats(:);
    end

    lons = lons_new;
    lats = lats_new;
end

axg = subplot(1,1,1);
hold on;
m_proj('miller', ...
       'long',[-180.0 180.0], ...
       'lat', [-87.5 87.5]);

m_pcolor(lons, lats, out_mat(:,:,1));
colorbar;
caxis([0, max(out_mat, [], 'all')]);

m_plotbndry('/usr/local/MATLAB/R2021a/toolbox/boundary/world', ...
    'color', 'k', 'linewidth', 1.5)
m_grid;
colormap(axg, flip(othercolor('Spectral6', 12)));

sgtitle("OMI/MLS climatology O3 column JANUARY - test [DU]");

set(gcf, 'Renderer', 'painters', 'Position', [90 0 1250 970]);
