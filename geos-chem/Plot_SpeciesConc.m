% Plot outputted file
% Haipeng Lin, 2020/06/08

opengl hardware;
colormap(flip(othercolor('Spectral10', 10)));

% for reusable code
data = transpose(tmp2D_O3(:,:,15)); % 09-25 00 = 1, 06 = 2, ...

% convert to ppb
% data = data * 1e9;

% below code
set(gcf, 'Renderer', 'painters', 'Position', [90 90 1000 600])

% [min(min(lons)) max(max(lons))]
% [min(min(lats)) max(max(lats))]
m_proj('miller', ...
       'long',[-135.0 -60.0], ...
       'lat', [10.0 80.0]);

m_pcolor(lons, lats, data);
hold on;
colorbar;
% caxis([1, 120]);

m_plotbndry('C:\Program Files\MATLAB\R2020a\toolbox\boundary\world', ...
    'color', 'k', 'linewidth', 1.5) % range lon is -180.0 180.0 so correct
m_plotbndry('C:\Program Files\MATLAB\R2020a\toolbox\boundary\shengjie', ...
    'color', 'k', 'linewidth', 1.5) % range lon is -180.0 180.0 so correct
m_grid;

t = title(strcat("Species O3 []"));
% Set font properties...
set(gca, 'FontName', 'Arial');
set(gca, 'FontWeight', 'normal');
set(gca, 'FontSize', 14);
set(findall(gcf,'type','text'), 'FontName', 'Arial');
set(findall(gcf,'type','text'), 'FontWeight', 'normal');
set(findall(gcf,'type','text'), 'FontSize', 14);

% Set printing properties
set(gcf, 'PaperUnits', 'inches', 'PaperPosition', [0 0 12 10]);