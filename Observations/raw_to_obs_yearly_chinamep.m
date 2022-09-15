% create yearly obs data from chinamep nc data

% first, build a dataset with all the observational time slices available
% i.e., 2016010100, 2016010101, ... 2016123124, 2017010108, ...
% consider that UTC+9 (KST) 2016010100(==2015123124) maps to 2015123115, etc., so this
% case needs to be handled.
%
% note that airkorea data has 24z instead of 00z; so data actually starts
% at 2016010101 which is 2015123116.

% (dates are awful to work with)
% can only process one year at a time.
year = 2016;     % edit here: year number
tz_offset = 8;   % edit here: 8 for UTC+8 (CST)

% data source format string - year, month
data_source_fstring = "raw_china_obs_nc_xfeng/%s_nonan_201301_202006.nc";
data_sink_fstring = "obs_chinamep_ground/%04d_hourly_all.mat";

% 24.45 or 22.4 (chinamep)
obs_gas_factor = 22.4;

% initialize data.
date_strings_list = [];
date_ptridx_list  = [];
curr_ptridx       = 1;

fprintf("Processing China MEP obs: tz_offset is UTC+%\n\n", tz_offset);

%% construct the final list of dates: date_strings_list and date2ptr
% fill periods surrounding tz_offset: part 1
if tz_offset > 0
    % need to rewind a few time slices; since no tz offset is more than 12h
    % it is always year-1, 12/31
    for i = 1:tz_offset
        % but invert this ... so we start filling from earliest

        curr_tz_offset = tz_offset - i + 1;
        time_tmp = 24-curr_tz_offset;
        date_string_tmp = sprintf("%04d1231%02d", year-1, time_tmp);
        date_strings_list = [date_strings_list, date_string_tmp];
        date_ptridx_list = [date_ptridx_list, curr_ptridx];
        curr_ptridx = curr_ptridx + 1;
    end
end

for month = 1:12
    rundays = 31;
    if month == 2
        rundays = 28;
        if (mod(year, 4) == 0 && mod(year, 100) ~= 0) || mod(year, 400) == 0
            rundays = 29;
        end
    elseif month == 4 || month == 6 || month == 9 || month == 11
        rundays = 30;
    end
    
    for day = 1:rundays
        for hour = 1:24
            date_string_tmp = sprintf("%04d%02d%02d%02d", year, ...
                                 month, day, hour);
            date_strings_list = [date_strings_list, date_string_tmp];
            date_ptridx_list = [date_ptridx_list, curr_ptridx];
            curr_ptridx = curr_ptridx + 1;
        end
    end
end

% fill periods surrounding tz_offset: part 2
% also advance a few time slices forward of year, to facilitate tz
% calculation and also to include all data avail
for i = 1:tz_offset
    time_tmp = i;
    date_string_tmp = sprintf("%04d0101%02d", year+1, time_tmp);
    date_strings_list = [date_strings_list, date_string_tmp];
    date_ptridx_list = [date_ptridx_list, curr_ptridx];
    curr_ptridx = curr_ptridx + 1;
end

date2ptr = containers.Map(date_strings_list, date_ptridx_list);
if tz_offset > 0
    final_num_idxs = curr_ptridx - 1 - tz_offset;
end

% for china mep
obs_vars_avail = ["no2", "so2", "o3", "pm25"];

%% read site data from nc (o3)
fname = sprintf(data_source_fstring, "o3");
% first build all the station data
siteids = ncread(fname, 'siteid');

% Resolve sites into site index, using another hash map
sitenum = length(siteids);
ids = 1:sitenum;
siteid2id = containers.Map(siteids, ids);

% sites: id: (siteid, lon, lat). retrieve using
% sites(id,:)=(siteid,lon,lat), it is a 3-dim array
sites = zeros(sitenum, 3);
sites(:,1) = siteids; % true site id
lons = ncread(fname, 'lon');
lats = ncread(fname, 'lat');
sites(:,2) = lons;
sites(:,3) = lats;

%% now go through each species and write time slices
for j = 1:length(obs_vars_avail)
    spcname = obs_vars_avail(j);
    fname = sprintf(data_source_fstring, spcname);

    current_fill_index = 1;
    if spcname == "pm25"
        data_raw = ncread(fname, 'pm2_5');
    else
        data_raw = ncread(fname, spcname);
    end

    % compute hour, day, month, year ptrs. note these are in CST,
    % to calculate this, simply lookahead z_offset timestamps forward
    for i = 1:final_num_idxs  % go through estimate time slices.
        datestring_cst = date_strings_list(i+tz_offset);
        datestring_utc = date_strings_list(i);
   
        % parse it out...
        cst_ptrs = sscanf(datestring_cst, "%04d%02d%02d%02d");

        % map the years to the appropriate num: -2012
        %     1 2 3 4 5 6 7 8
        %201  3 4 5 6 7 8 9 20
        cst_ptrs(1) = cst_ptrs(1) - 2012;

        % for china mep, convert ug/m3 to ppm for comparison
        convert_factor = 1;
        if spcname == "no2"
            convert_factor = obs_gas_factor * 1e-3 / 46.0055;
        elseif spcname == "o3"
            convert_factor = obs_gas_factor * 1e-3 / 48;
        elseif spcname == "so2"
            convert_factor = obs_gas_factor * 1e-3 / 64.066;
        end

        % note that obs_%s is first filled here, sites is not actually used
        eval(sprintf("obs_%s(:,current_fill_index) = data_raw(%d, %d, %d, %d, :) * convert_factor;", spcname, cst_ptrs(4), cst_ptrs(3), cst_ptrs(2), cst_ptrs(1)));
        fprintf("=> filled all stations at time %d/%d (cst=%s, utc=%s) species %s\n", current_fill_index, final_num_idxs, datestring_cst, datestring_utc, spcname);

        current_fill_index = current_fill_index + 1;
    end
end

% last filled index
fprintf("finished: last filled index is %d at %s (tz_offset=%d)\n\n", current_fill_index-1, ...
        date_strings_list(current_fill_index-1), tz_offset);

% post-process data with NaN
for j = 1:length(obs_vars_avail)
    spcname = obs_vars_avail(j);
    eval(sprintf("obs_%s(obs_%s < 0) = NaN;", spcname, spcname));
end

% save magic...
magic_string = "'date_strings_list', 'date2ptr', 'sitenum', 'sites', 'siteid2id'";

for j = 1:length(obs_vars_avail)
    eval(sprintf("obs_%s = obs_%s(:,1:%d);", obs_vars_avail(j), obs_vars_avail(j), current_fill_index-1));
    magic_string = append(magic_string, ", 'obs_", obs_vars_avail(j), "'");
end
magic_string

% last filled index
fprintf("finished: last filled index is %d at %s (tz_offset=%d)\n\n", current_fill_index-1, ...
        date_strings_list(current_fill_index-1), tz_offset);

magic_eval = sprintf("save('%s', %s)", sprintf(data_sink_fstring, year), magic_string);
magic_eval
eval(magic_eval);