function [tab] = fluxQuery(bucket, t0, tf, seriesname, tag_key, tag_value, field, resolution)
   % Define the IP, Port, and token
    IP = '128.179.34.35';
    PORT = 52000;
    token = "GqzP7J6acEsyH12LY9tf-be3r9FkOAgBWcL1ADCxFv8sP_D6r1kEMaARXsWjulepTAYrGfPHKjmV4zlC11PEPA==";
    org = "DESL-EPFL";
    url = sprintf('http://%s:%d/api/v2/query', IP, PORT);

    % Construct Flux query
    flux_query = sprintf(['from(bucket: "%s")\n' ...
                          '  |> range(start: %s, stop: %s)\n' ...
                          '  |> filter(fn: (r) => r["_measurement"] == "%s")\n' ...
                          '  |> filter(fn: (r) => r["%s"] == "%s")\n' ...
                          '  |> filter(fn: (r) => r["_field"] == "%s")\n' ...
                          '  |> aggregateWindow(every: %s,\n' ...
                          'fn: (tables=<-, column) => tables |> median(method: "exact_selector"))'], ...
                          bucket, t0, tf, seriesname, tag_key, tag_value, field, resolution);

    % Set headers with authorization token
    headers = {
        'Authorization', sprintf('Token %s', token);
        'Content-Type', 'application/json'
    };

    % Set query parameters
    queryParams = sprintf('org=%s&bucket=%s&precision=ms', org, bucket);

    % Convert the data structure to JSON
    data = struct('query', flux_query);
    jsonData = jsonencode(data);

    % Construct full URL with query parameters
    fullUrl = sprintf('%s?%s', url, queryParams);

    % Set web options for JSON request
    options = weboptions('Timeout',60,'HeaderFields', headers, 'MediaType', 'application/json', 'RequestMethod', 'post');

    % Make the HTTP request
    response = webwrite(fullUrl, jsonData, options);

    % Extract data from the table response
    if istable(response)
        % Assuming the response is a table
        time = datetime(response.x_time, 'InputFormat', 'yyyy-MM-dd''T''HH:mm:ss''Z''', 'TimeZone', 'UTC');
        data = response.x_value;
        
        tab = timetable(time,data,'VariableNames', field);
    else
        error('Unexpected response format.');
    end

end
