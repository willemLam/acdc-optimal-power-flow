
function [PV_perun, PV_solarmax, PV_facade] = PV_history(tf,t0)


database = 'logging';

seriesname = 'dayahead';
tag_keys = {'epfl_dayahead'};
tag_values = {'GHI'};

% [GHI_dayhead, ~] = readSeriesFromDatabase_InfluxDB2(database, seriesname, tag_keys, tag_values, [], t0, tf, 3600);

tag_values = {'airT_soda'};
[Tair_dayhead, ~] = readSeriesFromDatabase_InfluxDB2(database, seriesname, tag_keys, tag_values, [], t0, tf, 600);


seriesname = 'GHI';
seriesname = 'Global_Total_Irradiance';

tag_keys = {'device'};
tag_values = {'PVroof'};
[GHI_realization, ts] = readSeriesFromDatabase_InfluxDB2(database, seriesname, tag_keys, tag_values, [], t0, tf, 600);

% GHI_dayhead_last = nanmean(GHI_dayhead);
GHI_realization_last = nanmean(GHI_realization(end-2:end));

Tair_dayahead_last = nanmean(Tair_dayhead);
GHI_persistent = 0.9*(GHI_realization_last);
% - for persistent 
[PV_perun, PV_solarmax,PV_facade] = PV_5min_forecast_v3(GHI_persistent, Tair_dayahead_last, epochToMatlabTime(ts(end)), 2);
PV_perun(isnan(PV_perun)) = 0 ;
PV_solarmax(isnan(PV_solarmax)) = 0 ;
PV_facade(isnan(PV_facade)) = 0 ;

% save PV_forecast.mat PV_perun PV_solarmax
end

