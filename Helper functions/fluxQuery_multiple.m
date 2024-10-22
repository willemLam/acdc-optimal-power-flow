function [data, time] = fluxQuery(bucket, t0, tf, seriesname, tag_key, tag_value, field, resolution)
   % Define the IP, Port, and token
    IP = '128.179.34.35';
    PORT = 52000;
    token = "GqzP7J6acEsyH12LY9tf-be3r9FkOAgBWcL1ADCxFv8sP_D6r1kEMaARXsWjulepTAYrGfPHKjmV4zlC11PEPA==";
    org = "DESL-EPFL";
    url = sprintf('http://%s:%d/api/v2/query', IP, PORT);

    % Create the field filter for multiple fields
    fieldFilter = '';
    for i = 1:length(fields)
        if i == 1
            fieldFilter = sprintf('r["_field"] == "%s"', fields{i});
        else
            fieldFilter = sprintf('%s or r["_field"] == "%s"', fieldFilter, fields{i});
        end
    end

    % Construct Flux query with multiple fields
    flux_query = sprintf(['from(bucket: "%s")\n' ...
                          '  |> range(start: %s, stop: %s)\n' ...
                          '  |> filter(fn: (r) => r["_measurement"] == "%s")\n' ...
                          '  |> filter(fn: (r) => r["%s"] == "%s")\n' ...
                          '  |> filter(fn: (r) => %s)\n' ...
                          '  |> aggregateWindow(every: %s,\n' ...
                          'fn: (tables=<-, column) => tables |> median(method: "exact_selector"))'], ...
                          bucket, t0, tf, seriesname, tag_key, tag_value, fieldFilter, resolution);

    % Set headers with authorization token
    headers = {'Authorization', sprintf('Token %s', token);
               'Content-Type', 'application/json'};

    % Set query parameters
    params = sprintf('org=%s&bucket=%s&precision=s', org, bucket);
    fullUrl = sprintf('%s?%s', url, params);

    % Convert the data structure to JSON
    data = struct('query', flux_query);
    jsonData = jsonencode(data);
    
    % Set web options for JSON request
    options = weboptions('HeaderFields', headers, 'MediaType', 'application/json', 'RequestMethod', 'post');

    % Make the HTTP request
    response = webwrite(fullUrl, jsonData, options);

    % Extract data from the table response
    if istable(response)
        % Assuming the response is a table
        time = response.x_time;
        data = response.x_value;
    else
        error('Unexpected response format.');
    end

end
