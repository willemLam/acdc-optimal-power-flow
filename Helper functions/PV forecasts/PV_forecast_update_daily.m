function [PV_perun, PV_solarmax] = PV_forecast_update_daily()
UTC = 2;

tf = matlabToEpochTime(datenum(today+1))-UTC*3600;
t0 = matlabToEpochTime(datenum(today))-UTC*3600;

Hrs = hour(datestr(now));

tf = matlabToEpochTime(now)-UTC*3600 - Hrs*3600;
t0 = tf-24*3600;

database = 'logging';

seriesname = 'dayahead';
tag_keys = {'epfl_dayahead'};
tag_values = {'GHI'};

[GHI_dayhead, ts] = readSeriesFromDatabase_InfluxDB2(database, seriesname, tag_keys, tag_values, [], t0, tf, 3600);

tag_values = {'airT'};
[Tair_dayhead, ~] = readSeriesFromDatabase_InfluxDB2(database, seriesname, tag_keys, tag_values, [], t0, tf, 3600);

not_nan_idx = isnan(GHI_dayhead)~=1;
GHI_dayhead_notnan = GHI_dayhead(not_nan_idx);
Tair_dayahead_notnan = Tair_dayhead(not_nan_idx);
ts = ts(not_nan_idx);

[PV_perun, PV_solarmax] = PV_5min_forecast_v3(GHI_dayhead_notnan, Tair_dayahead_notnan, epochToMatlabTime(ts), 2);

PV_perun(isnan(PV_perun)) = 0 ;
PV_solarmax(isnan(PV_solarmax)) = 0 ;

% save PV_forecast.mat PV_perun PV_solarmax
end

