function [values, ts] = readSeriesFromDatabase_InfluxDB2( ...
    database, seriesname, tag_keys, tag_values, value,  t0, tf, resolution)
%READSERIESFROMDATABASE read <quantity> with the couple <tag_keys, tag_values>  from influxdb for the
% time period between <t0> and <tf> (epochtime, in second) at resolution
% <resolution> second.

% If the time series is empty or does not exist, it returns nan.
% to use this function:



if (nargin < 1)
%     % for debug only

    t0 = matlabToEpochTime(datenum(2024, 8, 7,19,1,48));
    tf = matlabToEpochTime(datenum(2024, 8, 7,19,5,59));
    seriesname = 'GHI';
    tag_keys = { 'device' };
    tag_values = { 'PVroof' };
    resolution = 1;
    database='microgrid_ST';


end

if isempty(value) == 1

    % Compose query
    groupbytime_clause = sprintf('time(%.0fs)', resolution);
    time_clause = sprintf('(time >= %.0fs AND time < %.0fs)', t0, tf);

    if numel(tag_keys) ~= numel(tag_values)
        error('numel(tag_keys) ~= numel(tag_values)!');
    end

    where_clause = sprintf('%s', time_clause);
    for i=1:numel(tag_keys)
        where_clause = sprintf('%s AND %s=''%s''', where_clause, tag_keys{i}, tag_values{i});
    end


    query = sprintf('SELECT mean(value) as value FROM %s WHERE %s GROUP BY %s', ...
        seriesname, where_clause, groupbytime_clause);

else
    
    % Compose query
    groupbytime_clause = sprintf('time(%.0fs)', resolution);
    time_clause = sprintf('(time >= %.0fs AND time < %.0fs)', t0, tf);

    if numel(tag_keys) ~= numel(tag_values)
        error('numel(tag_keys) ~= numel(tag_values)!');
    end

    where_clause = sprintf('%s', time_clause);
    for i=1:numel(tag_keys)
        where_clause = sprintf('%s AND %s=''%s''', where_clause, tag_keys{i}, tag_values{i});
    end


    query = sprintf(['SELECT mean(', value ') as value FROM %s WHERE %s GROUP BY %s'], ...
        seriesname, where_clause, groupbytime_clause);
    

end


% Perform query
% IP = '128.179.34.35';
IP = '128.179.34.35';
DB = "logging";

header = '--header "Authorization: Token GqzP7J6acEsyH12LY9tf-be3r9FkOAgBWcL1ADCxFv8sP_D6r1kEMaARXsWjulepTAYrGfPHKjmV4zlC11PEPA=="'; 
[exit_state, result] = system(sprintf('curl -X POST http://%s:52000/query?pretty=true -sS -H Accept:application/json > result_temp.json %s --data-urlencode db="%s" --data-urlencode q="%s"', IP, header, DB, query),'-echo') ;

% [a, b] = system(sprintf('curl -X POST "http://192.168.1.62:8086/query?pretty=true" -sS -H Accept:application/json > result.json --header "Authorization: Token X64jRFWeFhWFHHCkQvcDFPdmLVQX11l4aIamdObHNwnKFZW0YZQlTgTcBoswf1L5CeF-u0JJF5vzFmW2w_T7Mw==" --data-urlencode "db=logging" --data-urlencode "q=SELECT "GHI" FROM GITbox WHERE time > ''2021-02-01T10:19:59Z'' and time < ''2021-02-01T10:24:59Z'' AND "box"=''3''"'));


% html = urlread(url,'Timeout',3);
% html = strrep(html, 'null', 'nan');

result = loadjson('result_temp.json');


try
    n = max(size(result.results{1,1}.series{1, 1}.values));
    results = result.results{1,1}.series{1, 1}.values;
    values = results(:,2);
    ts = results(:,1);
catch
    values = nan;
    ts = nan;
end


if (nargin < 1)
    plot(ts, values)
    nop
end


if iscell(values)
    % Depending on the fact the json contains or not empty values,
    % the json has a different structure (shit!). This is a workaround to
    % still get the good values and setting to NaN the empty values.
    
%     values = [];
%     ts = [];

    clear ts values
    
    results = result.results{1,1}.series{1, 1}.values;
%     values = nan(n,1);
%     ts = nan(n,1);
    for i=1:n
        if ~iscell(results{i})
            if isnumeric(results{i}(1))
                ts(i) = results{i}(1);
                values(i) = results{i}(2);
            else
                ts_temp = (results{i}(1));
                ts_temp = char(strrep(ts_temp, 'T', ' '));
                ts(i,:) = ts_temp;
                values(i) = results{i}(2);
            end
        else
            ts_temp = (results{i}{1});
            ts_temp = char(strrep(ts_temp, 'T', ' '));
            ts(i,:) = ts_temp;
            values(i) = nan;
        end
    end
    
ts = matlabToEpochTime(datenum(ts(:,1:end-1)));

else
    
    
values = str2num(char(values));
ts = char(strrep(ts, 'T', ' '));

ts = datenum(ts(:,1:end-1));

    
end




end

