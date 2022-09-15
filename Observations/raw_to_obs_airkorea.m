% KORUS-AQ Project
% Process Korea Air Quality Monitoring Data, Ground
% Haipeng Lin <hplin@seas.harvard.edu>
% 2019-09-23 - Tell me if you wanna go home
% Updated 2022-09-12 - Numb/Encore - Linkin Park ft. Jay-Z
%
% Save data in "ground_tmp" should be a Nx9 table. 
% Use MATLAB - Home - Import Data to read raw CSV.
% Exclude row 1. Set datetime and site as strings. Import the rest as num.
% Replace blanks with NaN.
% Import as Table into "ground_tmp" var.
% 
%     Network      Site        Datetime       SO2     CO      O3      NO2     PM10    PM25
%     _______    ________    ____________    _____    ___    ____    _____    ____    ____
% 
%      City      "111121"    "2019050102"    0.003    0.5    0.02    0.044     21      11 

%% Read in the sites list from airkorea_stationinfo
% Size 347 x 3
sites = load("obs_airkorea_stationinfo.mat").AirKoreastationinfo1;
sitenum = size(sites,1);

% Resolve sites into site index, using another hash map
siteids = sites(:,1);
ids = 1:sitenum;
siteid2id = containers.Map(siteids, ids);

%% Generate timestamp hashmap
% ground_201905_temp(1,:).Datetime == "2019050101"
% timestamps go from 2019050101 to 2019050124 ... 2019053123
year = 2016;
month = 12;

% rundays can be dynamically set...
rundays = 31;
if month == 2
    rundays = 28;
    if (mod(year, 4) == 0 && mod(year, 100) ~= 0) || mod(year, 400) == 0
        rundays = 29;
    end
elseif month == 4 || month == 6 || month == 9 || month == 11
    rundays = 30;
end

ptr_index = 1;
strings_list = [];
ptrs_list = [];
for day = 1:rundays
    for hour = 1:24
        %if day == rundays && hour == 24
        %    break
        %end
        % also include 2019053124 because apparently it is part of May
        
        string_tmp = sprintf("%04d%02d%02d%02d", year, ...
                             month, day, hour);
        strings_list = [strings_list, string_tmp];
        ptrs_list = [ptrs_list, ptr_index];
        ptr_index = ptr_index + 1;
    end
end

date2ptr = containers.Map(strings_list, ptrs_list);
% test
% date2ptr("2016071224")

%% Now, read in the data row by row
% And store in array format
% sites: seqid: (siteid, lon, lat). retrieve using
% sites(seqid,:)=(siteid,lon,lat)
% obs: (seqid, timeptr)=spec_data. retrieve using obs_*(seqid,timeptr)

kor_obs_so2 = zeros(sitenum, ptr_index-1);
kor_obs_co = zeros(sitenum, ptr_index-1);
kor_obs_o3 = zeros(sitenum, ptr_index-1);
kor_obs_no2 = zeros(sitenum, ptr_index-1);
kor_obs_pm25 = zeros(sitenum, ptr_index-1);
kor_obs_pm10 = zeros(sitenum, ptr_index-1);

%f = waitbar(0,'Please wait...');
sz = size(ground_tmp, 1);

for i = 1:sz
    % if sidx does not exist, skip it
    if isKey(date2ptr, ground_tmp(i,:).Datetime)

        if class(ground_tmp(i,:).Site) == 'double'
            siteid = ground_tmp(i,:).Site;
        else
            siteid = str2num(ground_tmp(i,:).Site);
        end
    
        if isKey(siteid2id, siteid)
            sidx = siteid2id(siteid);
            didx = date2ptr(ground_tmp(i,:).Datetime);

            kor_obs_so2(sidx,didx) = ground_tmp(i,:).SO2;
            kor_obs_co(sidx,didx) = ground_tmp(i,:).CO;
            kor_obs_o3(sidx,didx) = ground_tmp(i,:).O3;
            kor_obs_no2(sidx,didx) = ground_tmp(i,:).NO2;
            kor_obs_pm25(sidx,didx) = ground_tmp(i,:).PM25;
            kor_obs_pm10(sidx,didx) = ground_tmp(i,:).PM10;
            if mod(i, 512) == 0
                fprintf("(%04d/%02d) processed %d out of %d observations\n", year, month, i, sz);
            end
        else
           %fprintf("Skipping STE %s out of range\n", ground_tmp(i,:).Site);
        end
        
    else
       %fprintf("Skipping %s out of range\n", ground_tmp(i,:).Datetime); 
    end
    
    %waitbar(i/sz,f,sprintf('%d/%d',i,sz))
end

%delete(f);

fprintf("Finished processing %04d/%02d total %d time slices\n", year, month, didx);

%% avoid namespace clash ... save to file
kor_sites = sites;
kor_siteid2id = siteid2id;
kor_sitenum = sitenum;
save(sprintf("obs_airkorea_ground/%04d_%02d_hourly.mat", year, month), "kor_siteid2id", "kor_sitenum", "kor_sites", ...
     "kor_obs_so2", "kor_obs_co", "kor_obs_o3", "kor_obs_no2", "kor_obs_pm25", "kor_obs_pm10");