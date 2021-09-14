% Plot O3 vertical profiles
% for CESM and GEOS-Chem Classic, hplin, 3/24/21

% ------ necessary input ------ %
% lons_GC, lats_GC (GC coords), EmisNO_Lightning_Col - GC var
% lons_CESM, lats_CESM (CESM coords), LNO_COL_PROD - CESM var
coords_GC = load('GC_GC_coords.mat');
lons_GC = coords_GC.lons;
lats_GC = coords_GC.lats;
levs_GC = coords_GC.levs .* 1000; % to hPa

coords_CESM = load('CESM_coords.mat');
lons_CESM = coords_CESM.lons;
lats_CESM = coords_CESM.lats;
levs_CESM = coords_CESM.levs;

clearvars coords_GC coords_CESM;

data_CESM = load('O3.mat');
O3_CESM = data_CESM.O3;

clearvars data_CESM;

data_GC = load('GC_O3.mat');
O3_GC = data_GC.GC_O3;

clearvars data_GC;

data_CAMC = load('CAMChem_O3.mat');
O3_CAMC = data_CAMC.O3;

clearvars data_CAMC;

%% plotting prs vs. CESM/GC
day_idx = 26;

z_range_GC = 1:72;
z_range_CESM = 1:56;
z_range_CESM_R = 56:-1:1;

% Find a station
% This assumes rectilinearity in xlat, xlong, also not very efficient as
% it goes to the end of the loop
slong = 127.310;
slat = 37.312;

% for CESM - note that CESM data is 0.0 to 360.0, and lons need to be 
% shifted to a -180.0 to 180.0 basis
if min(lons_CESM) < 0.0
    lons_CESM = lons_CESM - 180.0;
end

ddmax = 360;
di = 1;
dj = 1;
for i = 1:size(lons_CESM, 1)
    if abs(lons_CESM(i,1)-slong) < ddmax
        ddmsax = abs(lons_CESM(i,1)-slong);
        di = i;
    end
end

ddmax = 180;
for j = 1:size(lats_CESM, 1)
    if abs(lats_CESM(j,1)-slat) < ddmax
        ddmax = abs(lats_CESM(j,1)-slat);
        dj = j;
    end
end

i_CESM = di;
j_CESM = dj;
fprintf("i_CESM, j_CESM, %d, %d\n", i_CESM, j_CESM)

% for GC again
ddmax = 360;
di = 1;
dj = 1;
for i = 1:size(lons_GC, 1)
    if abs(lons_GC(i,1)-slong) < ddmax
        ddmsax = abs(lons_GC(i,1)-slong);
        di = i;
    end
end

ddmax = 180;
for j = 1:size(lats_GC, 1)
    if abs(lats_GC(j,1)-slat) < ddmax
        ddmax = abs(lats_GC(j,1)-slat);
        dj = j;
    end
end

i_GC = di;
j_GC = dj;
fprintf("i_GC, j_GC, %d, %d\n", i_GC, j_GC)

% mol/mol dry to ppbv
x_ozone_ppbv = reshape(O3_GC(i_GC, j_GC, z_range_GC, day_idx), [], 1) .* 1e9;
y_pressures = levs_GC(z_range_GC);

semilogy(x_ozone_ppbv, y_pressures, '-o');

ylabel("Pressure [hPa]");
xlabel("Ozone [ppb]");

hold on;

% plot CESM now
x_ozone_ppbv2 = reshape(O3_CESM(i_CESM, j_CESM, z_range_CESM, day_idx), [], 1) .* 1e9;
y_pressures2 = levs_CESM(z_range_CESM);

semilogy(x_ozone_ppbv2, y_pressures2, '-o');

% plot CAM-Chem now
x_ozone_ppbv3 = reshape(O3_CAMC(i_CESM, j_CESM, z_range_CESM, day_idx), [], 1) .* 1e9;
y_pressures3 = levs_CESM(z_range_CESM);

semilogy(x_ozone_ppbv3, y_pressures3, '-o');

% crude preprocessing; remove zeros...
x_obs_ozone_ppbv = [];
y_obs_prs = [];

for i = 1:size(Pressure_hPa, 1)/1.2
   if Ozone_ppmv(i) > 0
       x_obs_ozone_ppbv = [x_obs_ozone_ppbv; Ozone_ppmv(i) * 1000];
       y_obs_prs = [y_obs_prs; Pressure_hPa(i)];
   end
end

% plot obs sondes
%x_obs_ozone_ppbv = Ozone_ppmv .* 1000;
%y_obs_prs = Pressure_hPa;
plot(x_obs_ozone_ppbv, y_obs_prs, 'o', 'MarkerSize', 0.5);

% reverse y axis for pressure
set(gca, 'Ydir', 'reverse');

ylim([100 y_pressures(1)]);
%xlim([0 100]);

title(sprintf("Date: 2016-05-%02d, Lon: %f, Lat: %f", day_idx, slong, slat));

legend(["GEOS-Chem C", "CESM-GC", "CAM-Chem", "OBS (Taehwa Sonde)"], 'Location', 'southeast')

