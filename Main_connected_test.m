clear all
close all
addpath(genpath(pwd))
% addpath('C:\Users\admin\Documents\git-clones\PV forecasts\');
% addpath(genpath('C:\Users\admin\Documents\git-clones\PV forecasts'));
addpath(genpath('C:\Users\admin\Documents\MATLAB\YALMIP-master'));
addpath(genpath('C:\gurobi1002'));

%% RUN IN CMD
% cd C:\Users\admin\Documents\git-clones\acdc_centralized_controller2
% matlab -batch "Main_connected_test"

%% Initialization
tol = 1e-7; n_max = 100;
[Grid_para, Filter_para, idx1, idx3, constraints] = initialize_RT("CONNECTED");
idx = idx1; %single phase equivalent

% define control sampling time eg 30 sec
Ts = 0.5; %time in seconds

% Control parameters
TYPE = 'EXPERIMENT'; %'EXPERIMENT'; %'SIMULATION'; %
SEND_SETPOINTS = false;
SP_valid_termination = 0;
OPF_valid_termination = 0;

% Saving                          
FileName = fullfile(['C:\Users\admin\Documents\Willem\acdc_centralized_controller\results\',datestr(now,'mm-dd-yyyy HH-MM') ,'_connected_',TYPE,'.mat']);
emptyStruct = struct;
save(FileName,'-struct','emptyStruct');

% more initialization
N = 10000; % Length of experiment
yalmip('clear')
valid_termination = 0;
OPF_termination = 0;
prev_solution = 0;
% initialize for smooth start



%% Main control loop
for n_outer = 1:N-1

% Timing
T_outer = tic;
solution.time = datetime("now");

disp(['Control step ',num2str(n_outer)]);

% get data
variables = []; 
solution = Receive_raw_data_from_SCADA(Grid_para,variables,prev_solution);
[SM_info,variables] = ReceiveFromStateMachine_connected(Grid_para,variables) ;
disp(['Mode',SM_info.mode]);
disp([' ',[]]);

for n_inner = 1:4

    solution.Eabs = zeros(Grid_para.n_nodes,1);
    solution.P = zeros(Grid_para.n_nodes,1);
    solution.Q = zeros(Grid_para.n_nodes,1);

  

%     solution.Eabs(23:26) = constraints.Edc_nominal;%0.9375;
    solution.Eabs(23:26) = [constraints.Edc_nominal,0,0,constraints.Edc_nominal]; %constraints.Edc_nominal; %0.9375;
    solution.Eabs(27) = constraints.Edc_nominal;%0.9375;
    solution.Q(19:22) = [0,0.073,0.073,0.073];
    solution.P(27) = 0.0;
    solution.SupCap_P = 0;
    solution.Sambat_P = 0;
    solution.PV_P_perun = 0;
    solution.PV_P_solarmax = 0;

    prev_solution = solution;
    
end

    

%% Send the setpoints

constraints.Edc_nominal * Grid_para.Vdc_b
% To state machine
SM_info.timestep = posixtime(datetime("now"));
SM_info.feasibility = 1;%solution.proceed;
%AFE
setpoint_AFE = [solution.Eabs(23:26)*Grid_para.Vdc_b; solution.Q(19:22)*Grid_para.A_b]; 
% SUPERCAP
setpoint_SUPERCAP = [solution.Eabs(27)*Grid_para.Vdc_b; solution.SupCap_P*Grid_para.A_b];
%BESS
setpoint_BESS = [solution.Sambat_P*Grid_para.A_b/1000; 0];
%PERUN
setpoint_PERUN = [solution.PV_P_perun*Grid_para.A_b/1000; 0];
%SOLARMAX
setpoint_SOLARMAX = [solution.PV_P_solarmax*Grid_para.A_b/1000; 0];
%TOTAL
SM_info.setpoints = [setpoint_AFE;setpoint_SUPERCAP;setpoint_BESS;setpoint_PERUN;setpoint_SOLARMAX]';

SM_info = sendToStateMachine_connected(SM_info);

solution.duration = toc(T_outer);

%% Save all data
solution.duration = toc(T_outer);
disp(['Control finished in ', num2str(round(solution.duration,2)), ' seconds'])

name = ['data_',num2str(n_outer)];
assignin('base',name, solution)
save(FileName,name,'-append')
clear(genvarname(name))

disp('data saved');

disp('Time for a little nap!      zzz');

if TYPE == 'EXPERIMENT'
    time_to_wait = Ts - solution.duration;
    pause(time_to_wait);
%     
%     if OPF_termination >= 5
%         emergency_procedure(Grid_para)
%         error('The OPF didnt converge for more than 5 timesteps -> stop')
%     end
%     
%     if valid_termination >= 5
%         emergency_procedure(Grid_para)
%         error('The setpoint was not valid for more than 5 timesteps -> stop')
%     end
end

clear yalmip


end

