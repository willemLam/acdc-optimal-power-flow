function time_matlab = epochToMatlabTime(epoch_s, timezone)
% EPOCHTOMATLABTIME It converts from Unix epoch time (in seconds) to matlab time.
% TIMEZONE is the difference in hours between UTC and local time (e.g. +1 for CET)
% Fabrizio Sossan, DESL EPFL, 2016

% Matlab time is the number of days since 0 Jan 0000


    if nargin == 1
        timezone = 0;
    end
    
    epoch_s = epoch_s + 3600*timezone;
    offset = datenum('1970', 'yyyy');
    time_matlab = offset + epoch_s/8.64e4;
end
