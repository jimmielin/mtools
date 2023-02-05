% create netcdf for StateVector.nc
load('cluster_data_601_chn_separate.mat');
% Lat 160x1
% Lon 156x1
% clust_GMM 160x156 double

% create fileschema
ncSchema.Filename = 'StateVector.nc';
ncSchema.Name = '/';
ncSchema.Dimensions(1).Name = 'time';
ncSchema.Dimensions(1).Length = 1;
ncSchema.Dimensions(1).Unlimited = true;

ncSchema.Dimensions(2).Name = 'lat';
ncSchema.Dimensions(2).Length = size(Lat, 1);
ncSchema.Dimensions(2).Unlimited = false;

ncSchema.Dimensions(3).Name = 'lon';
ncSchema.Dimensions(3).Length = size(Lon, 1);
ncSchema.Dimensions(3).Unlimited = false;

ncSchema.Attributes(1).Name = 'Conventions';
ncSchema.Attributes(1).Value = 'CF-1.4';
ncSchema.Group = [];
ncSchema.Format = 'classic';
ncwriteschema('StateVector.nc', ncSchema);

% create the variables
latschema.Name                            = 'lat';
latschema.Dimensions.Name                 = 'lat';
latschema.Dimensions.Length               = length(Lat);
latschema.Dimensions.Unlimited            = false;
latschema.Datatype                        = 'single';
latschema.Format                          = 'classic';
latschema.Attributes(1).Name              = 'long_name';
latschema.Attributes(1).Value             = 'latitude';
latschema.Attributes(2).Name              = 'units';
latschema.Attributes(2).Value             = 'degrees_north';
ncwriteschema('StateVector.nc',latschema);
ncwrite('StateVector.nc','lat',Lat(:));

lonschema.Name                            = 'lon';
lonschema.Dimensions.Name                 = 'lon';
lonschema.Dimensions.Length               = length(Lon);
lonschema.Dimensions.Unlimited            = false;
lonschema.Datatype                        = 'single';
lonschema.Format                          = 'classic';
lonschema.Attributes(1).Name              = 'long_name';
lonschema.Attributes(1).Value             = 'longitude';
lonschema.Attributes(2).Name              = 'units';
lonschema.Attributes(2).Value             = 'degrees_east';
ncwriteschema('StateVector.nc',lonschema);
ncwrite('StateVector.nc','lon',Lon(:));

timeschema.Name                         = 'time';
timeschema.Dimensions.Name              = 'time';
timeschema.Dimensions.Length            = 1;
timeschema.Dimensions.Unlimited         = true;
timeschema.Datatype                     = 'double';
timeschema.Attributes(1).Name           = 'long_name';
timeschema.Attributes(1).Value          = 'simulation time';
timeschema.Attributes(2).Name           = 'units';
timeschema.Attributes(2).Value          = 'days since 2009-01-01 00:00:00';
timeschema.Attributes(3).Name           = 'calendar';
timeschema.Attributes(3).Value          = 'gregorian';
ncwriteschema('StateVector.nc',timeschema);
ncwrite('StateVector.nc','time',[0]);

% create variable...
varschema.Name                     = 'StateVector';
varschema.Dimensions(1).Name       = 'lon';
varschema.Dimensions(1).Length     = length(Lon);
varschema.Dimensions(1).Unlimited  = false;
varschema.Dimensions(2).Name       = 'lat';
varschema.Dimensions(2).Length     = length(Lat);
varschema.Dimensions(2).Unlimited  = false;
varschema.Dimensions(3).Name       = 'time';
varschema.Dimensions(3).Length     = 1; % only one slice
varschema.Dimensions(3).Unlimited  = true;
varschema.Datatype                 = 'single';
varschema.Attributes.Name          = 'units';
varschema.Attributes.Value         = '1';
ncwriteschema('StateVector.nc', varschema);
ncwrite('StateVector.nc', 'StateVector', clust_GMM');

ncinfo('StateVector.nc');