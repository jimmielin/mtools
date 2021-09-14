% Plot outputted file
% Haipeng Lin, 2020/06/08

opengl hardware;
colormap(jet(10));

% for reusable code
timeslice = 5+48+1; % n hours since (incl) 2016-07-01 00:00:00, total 93 hours until 20z
data1 = transpose(spc_1_NOx (:,:,timeslice)); % 09-25 00 = 1, 06 = 2, ...
data2 = transpose(spc_2x_NOx(:,:,timeslice));
% datadiff = (data2 - data1) ./ data1 * 100; % element-wise division for a % diff
datadiff = data2 - data1;

% convert to ppb
% data = data * 1e9;

% below code
set(gcf, 'Renderer', 'painters', 'Position', [90 90 1800 400])

% Set printing properties
set(gcf, 'PaperUnits', 'inches', 'PaperPosition', [0 0 18 4]);

%% create subplots
subplot(1, 3, 1); % mxn plot #

% [min(min(lons)) max(max(lons))]
% [min(min(lats)) max(max(lats))]
m_proj('miller', ...
       'long',[-155.0 -50.0], ...
       'lat', [10.0 80.0]);

m_pcolor(lons, lats, data1);
hold on;
colorbar;
caxis([0, 5e-9]);

m_plotbndry('C:\Program Files\MATLAB\R2020a\toolbox\boundary\world', ...
    'color', 'k', 'linewidth', 1.5) % range lon is -180.0 180.0 so correct
m_plotbndry('C:\Program Files\MATLAB\R2020a\toolbox\boundary\shengjie', ...
    'color', 'k', 'linewidth', 1.5) % range lon is -180.0 180.0 so correct
m_grid;

t = title(strcat("NO_x, 2 \times 2.5", char(176), " HEMCO [v/v dry]"));

% Set font properties...
set(gca, 'FontName', 'Arial');
set(gca, 'FontWeight', 'normal');
set(gca, 'FontSize', 14);
set(findall(gcf,'type','text'), 'FontName', 'Arial');
set(findall(gcf,'type','text'), 'FontWeight', 'normal');
set(findall(gcf,'type','text'), 'FontSize', 14);

%% subplot 2
subplot(1, 3, 2);
% [min(min(lons)) max(max(lons))]
% [min(min(lats)) max(max(lats))]
m_proj('miller', ...
       'long',[-155.0 -50.0], ...
       'lat', [10.0 80.0]);

m_pcolor(lons, lats, data2);
hold on;
colorbar;
caxis([0, 5e-9]);

m_plotbndry('C:\Program Files\MATLAB\R2020a\toolbox\boundary\world', ...
    'color', 'k', 'linewidth', 1.5) % range lon is -180.0 180.0 so correct
m_plotbndry('C:\Program Files\MATLAB\R2020a\toolbox\boundary\shengjie', ...
    'color', 'k', 'linewidth', 1.5) % range lon is -180.0 180.0 so correct
m_grid;

t = title(strcat("NO_x, 1 \times 1.25", char(176), " HEMCO [v/v dry]"));

% Set font properties...
set(gca, 'FontName', 'Arial');
set(gca, 'FontWeight', 'normal');
set(gca, 'FontSize', 14);
set(findall(gcf,'type','text'), 'FontName', 'Arial');
set(findall(gcf,'type','text'), 'FontWeight', 'normal');
set(findall(gcf,'type','text'), 'FontSize', 14);

%% subplot 2
subplot(1, 3, 3);
% [min(min(lons)) max(max(lons))]
% [min(min(lats)) max(max(lats))]
m_proj('miller', ...
       'long',[-155.0 -50.0], ...
       'lat', [10.0 80.0]);

m_pcolor(lons, lats, datadiff);
hold on;
colorbar;
caxis([-3.5e-10, 3.5e-10]);

m_plotbndry('C:\Program Files\MATLAB\R2020a\toolbox\boundary\world', ...
    'color', 'k', 'linewidth', 1.5) % range lon is -180.0 180.0 so correct
m_plotbndry('C:\Program Files\MATLAB\R2020a\toolbox\boundary\shengjie', ...
    'color', 'k', 'linewidth', 1.5) % range lon is -180.0 180.0 so correct
m_grid;

t = title(strcat("Exp - Ref"));

colormap(gca, (othercolor('BuDOr_12', 7)));

% Set font properties...
set(gca, 'FontName', 'Arial');
set(gca, 'FontWeight', 'normal');
set(gca, 'FontSize', 14);
set(findall(gcf,'type','text'), 'FontName', 'Arial');
set(findall(gcf,'type','text'), 'FontWeight', 'normal');
set(findall(gcf,'type','text'), 'FontSize', 14);

sgtitle('NO_x Concentration [2016-07-03 05z (01 EST, 22 PST)]')