%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%            _____            ______        %
% _______ _____  /_______________  /_______ %
% __  __ `__ \  __/  __ \  __ \_  /__  ___/ %
% _  / / / / / /_ / /_/ / /_/ /  / _(__  )  %
% /_/ /_/ /_/\__/ \____/\____//_/  /____/   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       "mtools" Research Toolkit           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Util_VarReader.m
%
% Reading arbitrary data slices from arbitrary model outputs.
% (c) 2019-2021 Haipeng Lin <jimmie.lin@gmail.com>
%
% Based off original WRF-GC script by xfeng <fengx7@pku.edu.cn>
%
% Version: 2021.09.28
% See changelog at end of script
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% mode: once|hourly|daily|index
mode           = 'index';

% if mode /= 'index', configure below:
%
% format: cesm|wrf|gcclassic
% select the format of date names used for output file names.
% CESM-style: yyyy-mm-dd-sssss
% WRF-style:  yyyy-mm-dd_HH:MM:SS
% GEOS-Chem Classic (gcclassic): GEOSChem.[CollectionName].yyyymmdd_HHII
%
% other formats can always be added.
%
% WARNING: An important limitation of this script is that
% one file must have EXACTLY one time slice THAT IS READ.
%
% For gcclassic, an extra feature is that "special" vars_in will
% automatically be allocated respective collection prefixes, i.e.
% Jval_O3O1D => GEOSChem.JValues, SpeciesConc_O3 => GEOSChem.SpeciesConc
% (TBD)
%
format         = 'cesm'; % wrf, cesm, gcclassic

% if mode /= 'index':
%
% filename: give an initial filename to read specs from,
% or the filename being read if mode == once
%
% if mode == 'index':
%
% filename: specify a FORMAT string where %d will be formatted
% to the index of given fslice_bounds. e.g., "test_%02d.nc", will be "test_01.nc",
% "test_02.nc", ...
filename       = 'FCSD_CAMChem.cam.h0.2016-2016._%02d_climo.nc';

% if mode == index: specify bounds for index (e.g. 1:12.) row vector only.
fslice_bounds  = 1:12;

% if mode /= 'once' && mode /= 'index':
%
% otherwise if mode /= 'once', provide the prefix and suffix
prefixname     = 'fcsd_camchem_clone.cam.h0.';
suffixname     = '.nc';

% if mode == 'hourly' || mode == 'daily', configure below date range:
rundays        = 30;
startyr        = 2016;
endyr          = 2016;
startmon       = 5;
endmon         = 5;
startdate      = 1; 
enddate        = 30;
num_days       = datenum(endyr,endmon,enddate) - datenum(startyr,startmon,startdate) + 1;

% fout_coords: output file name for coordinates spec
% will save lons, lats, levs
fout_coords = 'coords_CESM.mat';

% fout: output file name (all data will be saved to one file)
fout = 'climo_CAMchemYearlyMonthlyAvg_2016.mat';

% z_save [int]
% specify z-levels to save (for 3-D data), with level 1 being SURFACE
% ALL DATA WILL BE INVERTED EVEN IF TOA IS LEVEL 1, NO EXCEPTIONS
%
% if z_save > maximum available, then the maximum available will be saved.
% so if you want to save all, it is safe to set 999.
z_save = 56;

% variables to save
% vars_in_2D|3D will specify the variables that need to be read, with their
% original netCDF names...
%
% for many tools in mtools, the surface pressure and temperature are
% required. please store them in PSFC and T
vars_in_2D = ["PS", "LNO_COL_PROD", "DF_CO", "DF_O3"];
%vars_in_3D = ["O3", "NO", "NO2", "JvalO3O3P", "JvalO3O1D", "Jval_NO2", "Jval_CH2O"];
vars_in_3D = ["T", "O3", "NO", "NO2", "CO", "OH", "HO2", "SO2", "ISOP", "CH2O", "CLOUD", ...
    "BRO", "CLO", "HCL", "HOCL", "CH3CL", "N2O", "PAN"];%, "jo3_b", "jo3_a", "jno2", "jch2o_b", "jch2o_a"];

% vars_out_2D|3D will CORRESPONDINGLY rename these variables to be saved
% to the output (fout) file.
%
% BEWARE! ugly eval statements are used to save these variables.
% DO NOT USE ANY DANGEROUS MATLAB EXPRESSIONS HERE.
%
% (I have no business policing the variable names written below;
%  this is "on the other side of the airtight hatchway")
vars_out_2D = ["PSFC", "LNO_COL_PROD", "DryDepCO", "DryDepO3"];
vars_out_3D = ["T", "O3", "NO", "NO2", "CO", "OH", "HO2", "SO2", "ISOP", "CH2O", "CLDFRC" ...
    "BRO", "CLO", "HCL", "HOCL", "CH3CL", "N2O", "PAN"];% "JO3O3P", "JO3O1D", "JNO2", "JCH2O", "JCH2O_a"];

%%%%%%%%%%%%%%%%%%%%%% NO USER CONFIGURABLE CODE BELOW %%%%%%%%%%%%%%%%%%%%

% sanity check for input vars before we go further
if length(vars_in_2D) ~= length(vars_out_2D)
    if isempty(vars_out_2D) == false
        error("Error: size of vars_out_2D must be equal to vars_in_2D or be zero.");
    else
        fprintf("Warning: vars_out_2D is empty; setting it equal to vars_in_2D. You may want to use standardized names if comparing across models.\n");
    end
end

if length(vars_in_3D) ~= length(vars_out_3D)
    if isempty(vars_out_3D) == false
        error("Error: size of vars_out_3D must be equal to vars_in_3D or be zero.");
    else
        fprintf("Warning: vars_out_3D is empty; setting it equal to vars_in_3D. You may want to use standardized names if comparing across models.\n");
    end
end

% generate datestd list (from xfeng)
if strcmp(mode, 'daily')
    datestd        = zeros(num_days, 6);
elseif strcmp(mode, 'hourly')
    datestd        = zeros(24 * num_days, 6);
    % datestd dims: Y,M,D,H,I,S
elseif strcmp(mode, 'once')
    datestd  = zeros(1, 1);
    date_str = filename;
    fslices = 1;
elseif strcmp(mode, 'index')
    fslices = size(fslice_bounds, 2);
    filename_format_string = filename; % save this for later use
end
        
i = 0;
if strcmp(mode, 'daily') || strcmp(mode, 'hourly')
    for mm = datenum(startyr,startmon,startdate):datenum(endyr,endmon,enddate)
        if strcmp(mode, 'hourly')
            for hh = 0:23
                i = i + 1;
                datestd(i,:) = datevec(mm);
                datestd(i,4) = hh;
            end
        elseif strcmp(mode, 'daily')
            i = i + 1;
            datestd(i,:) = datevec(mm);
            datestd(i,4) = 0; % always 0-hour. change for slicing
        end
    end
    
    fslices = size(datestd,1);
end

% assumes 1 slice per file, 1 file per day
% (can also be 24 slices if using 'hourly')
if strcmp(format, 'cesm') == true
    if strcmp(mode, 'daily') == true
        formatOut = 'yyyy-mm-dd-00000';
        date_str = datestr(datestd,formatOut);
    end
    
    if strcmp(mode, 'hourly') == true
        % not implemented %
        error("format cesm and mode hourly is not supported (yet)");
    end
end

if strcmp(format, 'wrf') == true
    if strcmp(mode, 'once') == false
        formatOut = 'yyyy-mm-dd_HH:MM:SS';
        date_str = datestr(datestd,formatOut);
    end
end

%% read the coords
if strcmp(mode, 'index')
    filename_spec = sprintf(filename, fslice_bounds(1));
else
    filename_spec = filename;
end
lons = ncread(filename_spec, 'lon');
lats = ncread(filename_spec, 'lat');
levs = ncread(filename_spec, 'lev');

% check if we need to flip the vertical. this is so that TOA is always the
% highest indexed level, by convention in mtools
if levs(1) < levs(2)
    levs = flip(levs);
    do_flip_vert = true; 
    fprintf("Notice: levels were found as TOA = level 1 in input file. By convention in this set of tools, the vertical has been flipped so that coordinates face up.\n");
else
    do_flip_vert = false;
end

% check if units need to be corrected for the levs variable. in CESM,
% levs(1) is ~999 (101 kPa), unit is hPa. for GC, it is 1.000.
% TODO

% read hyai, hybi, hybm, hyam.
% warning: this may not always be available. TODO check
hyai = ncread(filename_spec, 'hyai');
hybi = ncread(filename_spec, 'hybi');
hyam = ncread(filename_spec, 'hyam');
hybm = ncread(filename_spec, 'hybm');

if do_flip_vert == true
    hyai = flip(hyai); hybi = flip(hybi);
    hyam = flip(hyam); hybm = flip(hybm);
end

if size(levs, 1) < z_save
    fprintf("WARNING: z_save > max levs. Setting as max available.\n")
    z_save = size(levs, 1);
end

save(fout_coords, "lons", "lats", "levs", "hyai", "hybi", "hyam", "hybm");
fprintf("Util_VarReader.m: Saved coordinates to %s\n", fout_coords);

%% read variables: 3-D vars
for i = 1:length(vars_in_3D)
    spiename = vars_in_3D(i);
    spiename_out = vars_out_3D(i);
    
    if strcmp(mode, 'index')
        tmp = ncread(filename_spec, spiename);
    else
        tmp = ncread(filename, spiename);
    end

    xdim = size(tmp,1);
    ydim = size(tmp,2);
    zdim = size(tmp,3);
    tmp = zeros(xdim,ydim,zdim,fslices);
    for tt = 1:fslices
        if strcmp(mode, 'daily') || strcmp(mode, 'hourly')
            filename = [prefixname,date_str(tt,:),suffixname];
        elseif strcmp(mode, 'index')
            filename = sprintf(filename_format_string, fslice_bounds(tt));
            fprintf("reading %s, fslice %d, %d\n", filename, tt, fslice_bounds(tt));
        end
        %disp(filename);
        tmp(:,:,:,tt) = ncread(filename,spiename);
    end
    
    if do_flip_vert == true
        tmp3D = tmp(:,:,zdim:-1:(zdim-z_save+1),:);
    else
        tmp3D = tmp(:,:,1:z_save,:);
    end
    
    % copy to out... DANGER
    eval(sprintf('%s=tmp3D;', spiename_out));

    fprintf("[ OK ] Read %s => %s, dim'l: %d x %d x %d\n", spiename, spiename_out, xdim, ydim, zdim);
end

for i = 1:length(vars_in_2D)
    spiename = vars_in_2D(i);
    spiename_out = vars_out_2D(i);
    
    if strcmp(mode, 'index')
        tmp = ncread(filename_spec, spiename);
    else
        tmp = ncread(filename, spiename);
    end

    xdim = size(tmp,1);
    ydim = size(tmp,2);
    tmp = zeros(xdim,ydim,fslices);
    for tt = 1:fslices
        if strcmp(mode, 'daily') || strcmp(mode, 'hourly')
            filename = [prefixname,date_str(tt,:),suffixname];
        elseif strcmp(mode, 'index')
            filename = sprintf(filename_format_string, fslice_bounds(tt));
        end
        %disp(filename);
        tmp(:,:,tt) = ncread(filename,spiename);
    end
    
    % copy to out... DANGER
    eval(sprintf('%s=tmp;', spiename_out));

    fprintf("[ OK ] Read %s => %s, dim'l: %d x %d (2D)\n", spiename, spiename_out, xdim, ydim);
end

if ~any(strcmp(vars_out_3D, 'T'))
    fprintf("\nWARNING: You did not save a temperature field 'T' in 3-D fields. If you use mtools to process data further, some computations may fail.\n")
    fprintf("Try to save out a temperature 3-D field in the 'T' variable.\n\n");
end

if ~any(strcmp(vars_out_2D, 'PSFC'))
    fprintf("\nWARNING: You did not save a surface pressure 'PSFC' field in 2-D fields. If you use mtools to process data further, some computations may fail.\n")
    fprintf("Try to save out a surface pressure 2-D field in the 'PSFC' variable.\n\n");
end

fprintf("Wait now... Writing files to disk.\n");

if isempty(vars_out_2D)
    save(fout, vars_out_3D{:});
elseif isempty(vars_out_3D)
    save(fout, vars_out_2D{:});
else
    save(fout, vars_out_2D{:}, vars_out_3D{:});
end

fprintf("Util_VarReader.m: Saved out fields to %s ... file slices %d\n", fout, fslices)
