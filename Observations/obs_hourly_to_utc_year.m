% concat hourly observational data in
% xxx_obs_yyy to yearly data with time slices

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
tz_offset = 9;   % edit here: 9 for UTC+9, etc.

% quirk: in airkorea data, data starts at 01 LT instead of 00 LT,
% which means that the next month's first hour of data is actually
% in the current slice. if this option is enabled, this quirk is
% accounted for.
obs_starts_at_01z = 1;

% data source format string - year, month
data_source_fstring = "obs_airkorea_ground/%04d_%02d_hourly.mat";
data_sink_fstring = "obs_airkorea_ground/%04d_hourly_all.mat";

% initialize data.
date_strings_list = [];
date_ptridx_list  = [];
curr_ptridx       = 1;

fprintf("Processing obs: tz_offset is UTC+%d \n\n", tz_offset);
if obs_starts_at_01z
    fprintf("*** QUIRK ENABLED: observations start at 01LT and end at 24LT in each monthly file\n")
end

%% construct the final list of dates: date_strings_list and date2ptr
% fill periods surrounding tz_offset: part 1
if tz_offset > 0
    % need to rewind a few time slices; since no tz offset is more than 12h
    % it is always year-1, 12/31
    for i = 1:tz_offset
        % but invert this ... so we start filling from earliest
        % quirk handling:
        %if i == 1 && obs_starts_at_01z
        %    continue
        %end

        time_tmp = 24 - tz_offset + i;
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

% actually just a guess
final_num_idxs = curr_ptridx - 1;

% peek into first file and construct final data arrays, as site list
% is unknown.

%% AirKorea: loop through monthly source data, and read in
for month = 1:12
    data_source = sprintf(data_source_fstring, year, month);

    if month == 1
        % peek into the first file and construct final data arrays
        % based on the first index size.
        data_source_meta = who('-file', data_source);
        if ismember("kor_sitenum", data_source_meta) % is airkorea data
            sitenum = load(data_source).kor_sitenum;
            sites = load(data_source).kor_sites;
            siteid2id = load(data_source).kor_siteid2id;
            import_mode = "airkorea";
        elseif ismember("chn_sitenum", data_source_meta) % is China MEP data
            % not implemented...
            import_mode = "chinamep";
        else
            import_mode = "unknown";
            data_source_meta
            error("Unrecognized data format in data_source.")
        end

        % initialize data structures. for airkorea, co/no2/o3/pm10/pm25/so2
        if import_mode == "airkorea"
            obs_vars_avail = ["co", "no2", "so2", "o3", "pm10", "pm25"];
            obs_vars_prefix = "kor";
        % china mep data: no2/o3/pm25/so2
        elseif import_mode == "chinamep"
            obs_vars_avail = ["no2", "so2", "o3", "pm25"];
            obs_vars_prefix = "chn";
        end

        % ugly ugly
        for j = 1:length(obs_vars_avail)
            eval(sprintf("obs_%s = zeros(sitenum, final_num_idxs);", obs_vars_avail(j)));
        end

        % start filling. assume data is continuous
        current_fill_index = 1;
    end

    % read the data from data_source based on import_mode.
    for j = 1:length(obs_vars_avail)
        eval(sprintf("tmp_obs_%s = load(data_source).%s_obs_%s;", obs_vars_avail(j), obs_vars_prefix, obs_vars_avail(j)));
    end

    tmp_obs_diml = size(tmp_obs_o3);
    time_slices_avail = tmp_obs_diml(2);

    % fill from current_fill_index, up to the maximum available time slices
    for j = 1:length(obs_vars_avail)
        eval(sprintf("obs_%s(:,%d:%d) = tmp_obs_%s;", obs_vars_avail(j), ...
            current_fill_index, current_fill_index+time_slices_avail-1, ...
            obs_vars_avail(j)));
        fprintf("%04d/%02d: filled %d time slices %d:%d for species %s\n", ...
            year, month, time_slices_avail, ...
            current_fill_index, current_fill_index+time_slices_avail-1, ...
            obs_vars_avail(j));
    end

    current_fill_index = current_fill_index + time_slices_avail;
end

% truncate the data because fills need to be absolutely correct
date_strings_list = date_strings_list(1:current_fill_index-1);
date_ptridx_list = date_ptridx_list(1:current_fill_index-1);

% fill metadata
date2ptr = containers.Map(date_strings_list, date_ptridx_list);
final_num_idxs = current_fill_index - 1;

% rename metadata (not needed, handled above)
% eval(sprintf("sitenum = %s_sitenum; sites = %s_sites; siteid2id = %s_siteid2id;", obs_vars_prefix, obs_vars_prefix, obs_vars_prefix));

% post-process data with NaN
for j = 1:length(obs_vars_avail)
    spcname = obs_vars_avail(j);
    eval(sprintf("obs_%s(obs_%s == 0) = NaN;", spcname, spcname));
end

% save magic...
magic_string = "'date_strings_list', 'date2ptr', 'sites', 'siteid2id', 'sitenum'";

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
