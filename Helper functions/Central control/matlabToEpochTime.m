function [ epoch_s, epoch_ms ] = matlabToEpochTime( matlab_time )
%MATLABTOEPOCHTIME Convert from matlab time to epoch time
% Fabrizio Sossan, DESL EPFL, 2016
  epoch_s = round(8.64e4 * (matlab_time - datenum(1970, 1, 1)));
  epoch_ms = round(8.64e7 * (matlab_time - datenum(1970, 1, 1)));
end