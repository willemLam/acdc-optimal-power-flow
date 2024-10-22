
function [PV_perun, PV_solarmax] = PV_forecast_update()
UTC = 2;
tf = matlabToEpochTime(now)-UTC*3600;
t0 = tf-2*3600;
database = 'logging';

seriesname = 'dayahead';
tag_keys = {'epfl_dayahead'};
tag_values = {'GHI'};

[GHI_dayhead, ~] = readSeriesFromDatabase_InfluxDB2(database, seriesname, tag_keys, tag_values, [], t0, tf, 3600);

tag_values = {'airT'};
[Tair_dayhead, ~] = readSeriesFromDatabase_InfluxDB2(database, seriesname, tag_keys, tag_values, [], t0, tf, 3600);


seriesname = 'GHI';
seriesname = 'Global_Total_Irradiance';

tag_keys = {'device'};
tag_values = {'PVroof'};
[GHI_realization, ts] = readSeriesFromDatabase_InfluxDB2_GHI(database, seriesname, tag_keys, tag_values, [], t0, tf, 60);

GHI_dayhead_last = nanmean(GHI_dayhead);
GHI_realization_last = nanmean(GHI_realization(end-2:end));

Tair_dayahead_last = nanmean(Tair_dayhead);
GHI_persistent = 0.9*(0*GHI_dayhead_last + 1*GHI_realization_last);
% - for persistent 
[PV_perun, PV_solarmax] = PV_5min_forecast_v3(GHI_persistent, Tair_dayahead_last, epochToMatlabTime(ts(end)), 2);
PV_perun(isnan(PV_perun)) = 0 ;
PV_solarmax(isnan(PV_solarmax)) = 0 ;

% save PV_forecast.mat PV_perun PV_solarmax
end

