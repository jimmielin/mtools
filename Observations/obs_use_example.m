% example to read observations

%% china mep

obs_chn_mep_2016 = load("obs_chinamep_ground/2016_hourly_all.mat");
obs_chn_mep_2016_sites = obs_chn_mep_2016.sites;
obs_chn_mep_2016_siteid2id = obs_chn_mep_2016.siteid2id;
obs_chn_mep_2016_date2ptr = obs_chn_mep_2016.date2ptr;
obs_chn_mep_2016_date_strings_list = obs_chn_mep_2016.date_strings_list;

obs_chn_mep_2016_o3 = obs_chn_mep_2016.obs_o3;

year = 2016; month = 2; day = 1; hour = 3; % utc time
seq_siteid = 1;
site = obs_chn_mep_2016_sites(seq_siteid, :);

seq_dateid = obs_chn_mep_2016_date2ptr(sprintf("%04d%02d%02d%02d", year, month, day, hour));

obs = obs_chn_mep_2016_o3(seq_siteid, seq_dateid);

fprintf("(mep) site: %d, lon: %f, lat: %f ", site(1), site(2), site(3));
fprintf("obs: o3 = %f\n", obs);

% sanity check against nc data
data_raw = ncread("raw_china_obs_nc_xfeng/o3_nonan_201301_202006.nc", "o3");
fprintf("raw obs: o3 = %f\n", data_raw(hour+8, day, month, year-2012, seq_siteid));


%% airkorea
obs_kor_airkorea_2016 = load("obs_airkorea_ground/2016_hourly_all.mat");
obs_kor_airkorea_2016_sites = obs_kor_airkorea_2016.sites;
obs_kor_airkorea_2016_siteid2id = obs_kor_airkorea_2016.siteid2id;
obs_kor_airkorea_2016_date2ptr = obs_kor_airkorea_2016.date2ptr;
obs_kor_airkorea_2016_date_strings_list = obs_kor_airkorea_2016.date_strings_list;

obs_kor_airkorea_2016_o3 = obs_kor_airkorea_2016.obs_o3;

year = 2016; month = 3; day = 1; hour = 4; % utc time
seq_siteid = 1;
site = obs_kor_airkorea_2016_sites(seq_siteid, :);

seq_dateid = obs_kor_airkorea_2016_date2ptr(sprintf("%04d%02d%02d%02d", year, month, day, hour));

obs = obs_kor_airkorea_2016_o3(seq_siteid, seq_dateid);

fprintf("(kor) site: %d, lon: %f, lat: %f ", site(1), site(2), site(3));
fprintf("obs: o3 = %f. check against raw data manually\n", obs);