function [ time_offset_h ] = getTimeZone(datestr)
%GETTIMEZONE Returns the timezone according to the local debian conf.
% It is the offset in hours between localtime and UTC/epoch.

t = datetime(datestr,'TimeZone','Europe/Zurich');

time_offset_h = hours(tzoffset(t));


end

