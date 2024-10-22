function [SM_info] = make_setpoint_array(SM_info,setpoint_AFE)

%% To state machine
SM_info.timestep = posixtime(datetime("now"));
SM_info.feasibility = 1;%solution.proceed;

%SUPERCAP
setpoint_SUPERCAP = [750; 0];
%BESS
setpoint_BESS = [0; 0];
%PERUN
setpoint_PERUN = [0; 0];
%SOLARMAX
setpoint_SOLARMAX = [0; 0];
%TOTAL
SM_info.setpoints = [setpoint_AFE;setpoint_SUPERCAP;setpoint_BESS;setpoint_PERUN;setpoint_SOLARMAX]';


end