function [ GNI_estimated ] = TransposedIrradiance(GHI, DN, UTC,panel_slope)

addpath(genpath('PV_LIB Version 1_32 Development'))
%GHI is the global horizontal irradiance 
%DN is date number obtained for example using datenum
%UTC is UTC time of Lausanne zone



%% Transposition model 
% panel_slope=10;
panel_azimuth=180;

%% Alternatively to grass we can use ineichen function 
%comment in case we want to use grass
Location.latitude = 46.518397;
Location.longitude = 6.565229;
Location.altitude = 200;

%Time
%DN = datenum(2017, 5 ,18, 0, 0,0):1/(24):datenum(2017, 5, 18, 23, 0, 0);

Time = pvl_maketimestruct(DN, UTC);


[ClearSkyGHI, ~, ClearSkyDHI]= pvl_clearsky_ineichen(Time, Location);

%we can get azimuth and elevation
[SunAz, SunEl, ~, ~] = pvl_ephemeris(Time, Location);


G_cs=ClearSkyGHI;
DHI=ClearSkyDHI;
G_b=G_cs-DHI;
AZIMUTH=SunAz;
ZENITH=90-SunEl;



for i=1:length(GHI)
if ZENITH(i)>80
GNI_estimated(i)=GHI(i);
else
GNI_estimated(i) = getTransposedIrradiance( ...
GHI(i), G_cs(i), G_b(i),  panel_slope, panel_azimuth,  ZENITH(i), AZIMUTH(i));
end
end

