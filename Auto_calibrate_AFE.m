clear all
close all
addpath(genpath(pwd))
% addpath('C:\Users\admin\Documents\git-clones\PV forecasts\');
% addpath(genpath('C:\Users\admin\Documents\git-clones\PV forecasts'));
addpath(genpath('C:\Users\admin\Documents\MATLAB\YALMIP-master'));
addpath(genpath('C:\gurobi1002'));
warning('off','all')

%% RUN IN CMD
% cd C:\Users\admin\Documents\git-clones\acdc_centralized_controller_3AFE
% matlab -batch "Main"

%% Initialization
tol = 1e-7; n_max = 100;
[Grid_para, Filter_para, idx1, idx3, constraints] = initialize_RT("CONNECTED");
idx = idx1;
variables = [];
variables = Receive_data_from_SCADA(Grid_para,variables,constraints,idx,"CONNECTED");

SM_info = [];
[SM_info,variables] = ReceiveFromStateMachine(SM_info,Grid_para,variables,"CONNECTED");

idx = idx1; %single phase equivalent


%% Get the data
start_time = tic();
ranges = [-7000:1000:7000];
n = 3;

% data.IC_1_Q = zeros(length(ranges));
% data.IC_2_Q = zeros(length(ranges));
% data.IC_3_Q = zeros(length(ranges));
% 
% data.IC_2_P = zeros(length(ranges));
% data.IC_3_P = zeros(length(ranges));

%% IC_2_P
for i = 1:length(ranges)
    ranges(i)
    setpoint_AFE = [constraints.Edc_nominal*Grid_para.Vdc_b; ranges(i) ; 0; constraints.Edc_nominal*Grid_para.Vdc_b; 0;0;0 ; 0]; 
    SM_info = make_setpoint_array(SM_info,setpoint_AFE);
    sendToStateMachine(SM_info,"CONNECTED");

    pause(n)

    variables = Receive_data_from_SCADA(Grid_para,variables,constraints,idx,"CONNECTED");
    data_IC_2_P(i) = real(variables.S(20))*Grid_para.A_b;
    sp_IC_2_P(i) = ranges(i);
end
setpoint_AFE = [constraints.Edc_nominal*Grid_para.Vdc_b; 0 ; 0; constraints.Edc_nominal*Grid_para.Vdc_b; 0;0;0 ; 0]; 
SM_info = make_setpoint_array(SM_info,setpoint_AFE);
sendToStateMachine(SM_info,"CONNECTED");

fit_data = fitlm(data_IC_2_P, sp_IC_2_P ).Coefficients;

figure
scatter(ranges ,data_IC_2_P )
title(['IC 2 - P - linear fit: ', num2str(fit_data.Estimate(2)),' x + ', num2str(fit_data.Estimate(1))])

%% IC_3_P
for i = 1:length(ranges)
    ranges(i)
    setpoint_AFE = [constraints.Edc_nominal*Grid_para.Vdc_b; 0 ; ranges(i) ; constraints.Edc_nominal*Grid_para.Vdc_b; 0;0;0 ; 0]; 
    SM_info = make_setpoint_array(SM_info,setpoint_AFE);
    sendToStateMachine(SM_info,"CONNECTED");

    pause(n)

    variables = Receive_data_from_SCADA(Grid_para,variables,constraints,idx,"CONNECTED");
    data_IC_3_P(i) = real(variables.S(21))*Grid_para.A_b;
    sp_IC_3_P(i) = ranges(i);
end
setpoint_AFE = [constraints.Edc_nominal*Grid_para.Vdc_b; 0 ; 0; constraints.Edc_nominal*Grid_para.Vdc_b; 0;0;0 ; 0]; 
SM_info = make_setpoint_array(SM_info,setpoint_AFE);
sendToStateMachine(SM_info,"CONNECTED");

fit_data = fitlm(data_IC_3_P ,sp_IC_3_P).Coefficients;

figure
scatter(ranges ,data_IC_3_P )
title(['IC 3 - P - linear fit: ', num2str(fit_data.Estimate(2)),' x + ', num2str(fit_data.Estimate(1))])


%% IC_1_Q
for i = 1:length(ranges)
    ranges(i)
    setpoint_AFE = [constraints.Edc_nominal*Grid_para.Vdc_b; 0 ; 0 ; constraints.Edc_nominal*Grid_para.Vdc_b; ranges(i); 0;0 ; 0]; 
    SM_info = make_setpoint_array(SM_info,setpoint_AFE);
    sendToStateMachine(SM_info,"CONNECTED");

    pause(n)

    variables = Receive_data_from_SCADA(Grid_para,variables,constraints,idx,"CONNECTED");
    data_IC_1_Q(i) = imag(variables.S(19))*Grid_para.A_b;
    sp_IC_1_Q(i) = ranges(i);
end
setpoint_AFE = [constraints.Edc_nominal*Grid_para.Vdc_b; 0 ; 0; constraints.Edc_nominal*Grid_para.Vdc_b; 0;0;0 ; 0]; 
SM_info = make_setpoint_array(SM_info,setpoint_AFE);
sendToStateMachine(SM_info,"CONNECTED");

fit_data = fitlm(data_IC_1_Q , sp_IC_1_Q).Coefficients;

figure
scatter(ranges ,data_IC_1_Q )
title(['IC 1 - Q - linear fit: ', num2str(-fit_data.Estimate(2)),' x + ', num2str(fit_data.Estimate(1))])


%% IC_2_Q
for i = 1:length(ranges)
    ranges(i)
    setpoint_AFE = [constraints.Edc_nominal*Grid_para.Vdc_b; 0 ; 0; constraints.Edc_nominal*Grid_para.Vdc_b; 0; ranges(i) ;0 ; 0]; 
    SM_info = make_setpoint_array(SM_info,setpoint_AFE);
    sendToStateMachine(SM_info,"CONNECTED");

    pause(n)

    variables = Receive_data_from_SCADA(Grid_para,variables,constraints,idx,"CONNECTED");
    data_IC_2_Q(i) = imag(variables.S(20))*Grid_para.A_b;
    sp_IC_2_Q(i) = ranges(i);
end
setpoint_AFE = [constraints.Edc_nominal*Grid_para.Vdc_b; 0 ; 0; constraints.Edc_nominal*Grid_para.Vdc_b; 0;0;0 ; 0]; 
SM_info = make_setpoint_array(SM_info,setpoint_AFE);
sendToStateMachine(SM_info,"CONNECTED");

fit_data = fitlm(data_IC_2_Q,sp_IC_2_Q).Coefficients

figure
scatter(ranges,data_IC_2_Q)
title(['IC 2 - Q - linear fit: ', num2str(-fit_data.Estimate(2)),' x + ', num2str(fit_data.Estimate(1))])


%% IC_3_Q
for i = 1:length(ranges)
    ranges(i)
    setpoint_AFE = [constraints.Edc_nominal*Grid_para.Vdc_b; 0; 0; constraints.Edc_nominal*Grid_para.Vdc_b; 0; 0; ranges(i); 0]; 
    SM_info = make_setpoint_array(SM_info,setpoint_AFE);
    sendToStateMachine(SM_info,"CONNECTED");

    pause(n)

    variables = Receive_data_from_SCADA(Grid_para,variables,constraints,idx,"CONNECTED");
    data_IC_3_Q(i) = imag(variables.S(21))*Grid_para.A_b;
    sp_IC_3_Q(i) = ranges(i);
end
setpoint_AFE = [constraints.Edc_nominal*Grid_para.Vdc_b; 0 ; 0; constraints.Edc_nominal*Grid_para.Vdc_b; 0;0;0 ; 0]; 
SM_info = make_setpoint_array(SM_info,setpoint_AFE);
sendToStateMachine(SM_info,"CONNECTED");

fit_data = fitlm(data_IC_3_Q,sp_IC_3_Q).Coefficients;

figure
scatter(sp_IC_3_Q,data_IC_3_Q)
title(['IC 3 - Q - linear fit: ', num2str(-fit_data.Estimate(2)),' x + ', num2str(fit_data.Estimate(1))])


disp(['the calibration took ', num2str(toc(start_time)),' seconds'])


%% process the data




