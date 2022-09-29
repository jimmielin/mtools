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
% Also plots O3 column against OMI/MLS climatology
% if the O3 species is available.
%
% Also can read a set of reference mass values from
% Util_ReadReferenceGlobalMass.m if available,
% used to compare against GEOS-Chem benchmarking
% results in FullChem/BenchmarkResults/Tables/\
% GlobalMass_{Trop,TropStrat}_{Jan,Apr,Jul,Oct}2019.txt
%
% (c) 2019-2022 Haipeng Lin <jimmie.lin@gmail.com>
%
% Version: 2022.09.15
% Started: 2022.09.15
% See changelog at end of script
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

data_source = "raw_gc_benchmark/13.4.0-rc.4-GlobalMass_TropStrat_Jul2019.txt";
data_sink = "mdl_ref_gc_benchmark/mdl_ref_gc_benchmark_13.4.0rc.4-tropstrat-201907.mat";

% create the target map
mdl_ref_sum_map = containers.Map;

% fill one by one
fid = fopen(data_source);

tline = fgetl(fid);
while ischar(tline)
    % advance one line... (note this means first line will be skipped)
    tline = fgetl(fid);

    if tline == -1
        break
    end

    % match using magic
    [tokens, matches] = regexp(tline, "^(\w+)\s*:\s*[^\s]+\s*(\d+\.\d{6})", "tokens", "match");
    if size(tokens, 2) == 0
        continue
    end

    mdl_ref_sum_map(tokens{1}{1}) = str2double(tokens{1}{2});
    fprintf("=> %10s is %4.6f Gg\n", tokens{1}{1}, str2double(tokens{1}{2}));
end

fclose(fid);

% save out
save(data_sink, "mdl_ref_sum_map");
fprintf("completed import of gc benchmark mass, sink = %s...\n", data_sink);