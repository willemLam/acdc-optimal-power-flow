function [data, time] = fluxQuery(bucket, t0, tf, seriesname, tag_key, tag_value, field, resolution)
    % fluxQuery sends a Flux query to an InfluxDB instance and returns the data and time values.
    % 
    % Inputs:
    %   bucket      - The InfluxDB bucket name
    %   t0          - Start time in RFC3339 format (e.g., '2021-02-01T10:00:00Z')
    %   tf          - End time in RFC3339 format (e.g., '2021-02-01T11:00:00Z')
    %   seriesname  - The measurement name
    %   tag_key     - The tag key to filter by
    %   tag_value   - The tag value to filter by
    %   field       - The field key to filter by
    %   resolution  - The aggregation resolution (e.g., '1m', '5m', '1h')
    %
    % Outputs:
    %   data        - The extracted data values
    %   time        - The corresponding time values
    
    IP = '128.179.34.35';
    PORT = 52000;
    token = "qfMF0UgjLmS3iQUtwFm6nZPIkRtpp-Xz2PZGBsWMG_jxkZhay2bIfIOd6l3Cd2QB0Pd0KQkws9bXXR1FFIK16Q==";
    org = "DESL-EPFL";
    url = sprintf('http://%s:%d/api/v2/query', IP, PORT);

    % Construct the Flux query
    flux_query = sprintf([...
        'from(bucket: "%s") |> range(start: %s, stop: %s) ', ...
        '|> filter(fn: (r) => r["_measurement"] == "%s") ', ...
        '|> filter(fn: (r) => r["%s"] == "%s") ', ...
        '|> filter(fn: (r) => r["_field"] == "%s") ', ...
        '|> aggregateWindow(every: %s, fn: (tables=<-, column) => tables |> median(method: "exact_selector"))'], ...
        bucket, t0, tf, seriesname, tag_key, tag_value, field, resolution);

    % Set the headers
    headers = {
        '-H', sprintf('Authorization: Token %s', token), ...
        '-H', 'Content-Type: application/json'
    };

    % Set the query parameters and data payload
    params = {
        '-d', sprintf('{"query": "%s", "type": "flux", "org": "%s", "bucket": "%s", "precision": "s"}', ...
                      flux_query, org, bucket)
    };

    % Execute the curl command
    curl_cmd = sprintf('curl -X POST %s %s %s', strjoin(headers), url, strjoin(params));
    [status, cmdout] = system(curl_cmd);

    % Error handling
    if status ~= 0
        error('Failed to execute Flux query: %s', cmdout);
    end

    % Parse the response
    lines = strsplit(cmdout, '\n');
    lines = lines(~cellfun('isempty', lines)); % Remove empty lines
    n = numel(lines);

    % Initialize output arrays
    data = zeros(n-1, 1);
    time = strings(n-1, 1);

    for i = 2:n
        fields = strsplit(lines{i}, ',');
        time(i-1) = fields{6}; % Adjust index based on actual response structure
        data(i-1) = str2double(fields{7}); % Adjust index based on actual response structure
    end

    % Convert time strings to datetime
    time = datetime(time, 'InputFormat', 'yyyy-MM-dd''T''HH:mm:ss.SSSSSS''Z', 'TimeZone', 'UTC');

end


