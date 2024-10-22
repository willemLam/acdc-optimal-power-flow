clear all
close all
addpath(genpath(pwd))
% addpath('C:\Users\admin\Documents\git-clones\PV forecasts\');
% addpath(genpath('C:\Users\admin\Documents\git-clones\PV forecasts'));
addpath(genpath('C:\Users\admin\Documents\MATLAB\YALMIP-master'));
addpath(genpath('C:\gurobi1002'));
warning('off','all')

%% RUN IN CMD
% cd C:\Users\admin\Documents\git-clones\acdc_centralized_controller_3AFE_PQ
% matlab -batch "Main_island"

%% Initialization
tol = 1e-7; n_max = 20;
[Grid_para, Filter_para, idx1, idx3, constraints] = initialize_RT("ISLAND");
idx = idx1; %single phase equivalent


% define control sampling time eg 30 sec
Ts = 1; %time in seconds

% Control parameters
TYPE = 'EXPERIMENT'; %'EXPERIMENT'; %'SIMULATION'; %
LIMIT_sp = false;
SEND_SETPOINTS = false;
SP_valid_termination = 0;
OPF_valid_termination = 0;

% Saving                          
FileName = fullfile(['C:\Users\admin\Documents\Willem\acdc_centralized_controller\results\',datestr(now,'mm-dd-yyyy HH-MM') ,'_',TYPE,'_island.mat']);
emptyStruct = struct;
save(FileName,'-struct','emptyStruct');

% more initialization
N = 10000; % Length of experiment
yalmip('clear')
valid_termination = 0;
OPF_termination = 0;

% initialize for smooth start
variables = [];
prev_solution = [];
SM_info = [];

variables = Receive_data_from_SCADA(Grid_para,variables,constraints,idx,"ISLAND");



prev_solution.P = variables.P;
prev_solution.Q = variables.Q;
prev_solution.S = variables.S;
prev_solution.Eabs = variables.Eabs;
prev_solution.Eang = variables.Eang;
prev_solution.E = prev_solution.Eabs.*exp(1i*prev_solution.Eang);
prev_solution.Iabs = variables.Iabs;
prev_solution.I = variables.I;


%% Main control loop
for n_outer = 1:N-1

% Timing
T_outer = tic;
time_loop = datetime("now")

disp(['Control step ',num2str(n_outer)]);

%% get data
variables = Receive_data_from_SCADA(Grid_para,variables,constraints,idx,"ISLAND");

%% offline
if TYPE == 'SIMULATION'
    if n_outer > 1
        solution.reconstructed.S(3) = -0.12  ;
        solution.reconstructed.S(11) = 0.15  ;
        tol = 1e-7; n_max = 100;
        [solution.reconstructed.E,~,~,~] = NR_rectangularACDC_1ph_general_V2_quadratic_loss(Grid_para,Filter_para,solution.reconstructed.S,solution.reconstructed.E,solution.reconstructed.E,idx,tol,n_max);
        solution.reconstructed.S = solution.reconstructed.E.*conj(Grid_para.YY*solution.reconstructed.E);
        solution.reconstructed.I = get_Current_flow(solution.reconstructed.E,Grid_para);
        variables.E = solution.reconstructed.E;
        variables.Eabs = abs(variables.E);
        variables.Eang = angle(variables.E);
        variables.S = solution.reconstructed.S;
        variables.P = real(variables.S);
        variables.Q = imag(variables.S);
        variables.I = get_Current_flow(solution.reconstructed.E,Grid_para);
        variables.Iabs = abs(variables.I);
    end
end

[SM_info,variables] = ReceiveFromStateMachine(SM_info,Grid_para,variables,"ISLAND");
disp(['Mode',SM_info.mode]);
disp([' ',[]]);

% variables.S(26) =  0.01;
% variables.S(20) =  variables.S(20) + 0.02;
% variables.S(21) =  variables.S(21) - 0.03;
% [E,~,~] = NR_rectangularACDC_1ph_general_V2_quadratic_loss(Grid_para,Filter_para,variables.S,variables.E,variables.E,idx1,tol,n_max);
% S = E.*conj(Grid_para.YY*E);
% Iac = abs(Grid_para.YY(Grid_para.pos_ac3(:,1),:) * E);
% Idc = abs(Grid_para.YY(Grid_para.pos_dc3(:,1),:) * E);
% Vdc = abs(E(Grid_para.pos_dc3(:,1)));
% P_quad_ic = Filter_para.a .* Vdc + Filter_para.b .* Iac + Filter_para.c .* Iac.^2 ...
%                                  + Filter_para.d .* Idc + Filter_para.e .* Idc.^2

%%

for n_inner = 1:3

    if n_inner == 3 && toc(T_outer) > 0.75
        continue
    end

    if n_inner == 1
        temp_variables = variables;
        state = variables;
    else
        temp_variables.E = solution.reconstructed.E;
        temp_variables.S = solution.reconstructed.S;
        temp_variables.Iabs = solution.reconstructed.Iabs;
        temp_variables.P = real(temp_variables.S);
        temp_variables.Q = imag(temp_variables.S);
        temp_variables.Eabs = abs(temp_variables.E);
        temp_variables.Eang = angle(temp_variables.E);
                
    end

    % OPF
    step = exp(1-n_inner);
    [solution,OPF_termination] = make_optimizer_island_V3_3(step,SM_info,Grid_para,Filter_para,idx,constraints,temp_variables,prev_solution,OPF_termination,Ts,state);
    % load flow
    [solution.reconstructed.E,~,~] = NR_rectangularACDC_1ph_general_V2_quadratic_loss(Grid_para,Filter_para,solution.S,solution.E,solution.E,idx1,tol,n_max);
    solution.reconstructed.S = solution.reconstructed.E.*conj(Grid_para.YY*solution.reconstructed.E);
    solution.reconstructed.Iabs = abs(get_Current_flow(solution.reconstructed.E,Grid_para));

    fprintf('Difference OPF and LF    %f p.u. \n', max(abs(abs(solution.reconstructed.E) - solution.Eabs)))
    
end

%  solution.S2 = solution.E.*conj(Grid_para.YY*solution.E);
%  Iac = abs(Grid_para.YY(Grid_para.pos_ac3(:,1),:) * solution.E);
%  Idc = abs(Grid_para.YY(Grid_para.pos_dc3(:,1),:) * solution.E);
%  Vdc = abs(solution.E(Grid_para.pos_dc3(:,1)));
%  P_quad_ic = Filter_para.a .* Vdc + Filter_para.b .* Iac + Filter_para.c .* Iac.^2 ...
%                                   + Filter_para.d .* Idc + Filter_para.e .* Idc.^2



%% limit the setpoint variations
if n_outer ~= 1 && LIMIT_sp == true
    V_limit = 5/Grid_para.Vdc_b; %5V 
    P_limit = 1000/Grid_para.A_b; %1000W

    solution.PV_P_solarmax = setpoint_limit(solution.PV_P_solarmax,prev_solution.PV_P_solarmax,P_limit);
    solution.PV_P_perun = setpoint_limit(solution.PV_P_perun,prev_solution.PV_P_perun,P_limit);
    solution.Sambat_P = setpoint_limit(solution.Sambat_P,prev_solution.Sambat_P,2*P_limit);
    solution.SupCap_P = setpoint_limit(solution.SupCap_P,prev_solution.SupCap_P,P_limit);
    solution.IC_Q = setpoint_limit(solution.IC_Q,prev_solution.IC_Q,P_limit);
    solution.IC_E = setpoint_limit(solution.IC_E,prev_solution.IC_E,V_limit);

    solution.E(24:26) = solution.IC_E;
    solution.S(19:22) = complex(solution.P(19:22),solution.IC_Q);
    solution.S(9) = complex(solution.PV_P_perun + solution.Sambat_P,solution.Q(9));
    solution.S(11) = complex(solution.PV_P_solarmax,solution.Q(11));
    solution.S(27) = solution.SupCap_P;

    [solution.reconstructed.E,~,~] = NR_rectangularACDC_1ph_general_V2(Grid_para,Filter_para,solution.S,solution.E,solution.E,idx1,tol,n_max);
    solution.reconstructed.S = solution.reconstructed.E.*conj(Grid_para.YY*solution.reconstructed.E);
    solution.reconstructed.Iabs = abs(get_Current_flow(solution.reconstructed.E,Grid_para));
end


%% Print the results
Show_setpoints(solution,variables,constraints,Grid_para,idx);


%% Validate the setpoints
solution.proceed = CheckSetpoints(solution, constraints, Grid_para, idx);

if solution.proceed == 1
    valid_termination = 0;
    prev_solution = solution;

elseif solution.proceed == 0
    warning('At least one setpoint is not valid -> Previous setpoint is implemented')
    solution = prev_solution;
    solution.proceed = 0;

    valid_termination = valid_termination+1;
end 
    
%% Send the setpoints

% To state machine
SM_info.timestep = posixtime(datetime("now"));
SM_info.feasibility = 1;%solution.proceed;
%AFE
setpoint_AFE = [constraints.Edc_nominal*Grid_para.Vdc_b; solution.IC_P*Grid_para.A_b; constraints.Edc_nominal*Grid_para.Vdc_b; -solution.IC_Q*Grid_para.A_b; 0];
%SUPERCAP
setpoint_SUPERCAP = [solution.Eabs(27)*Grid_para.Vdc_b; -solution.SupCap_P*Grid_para.A_b]
%BESS
setpoint_BESS = [solution.Sambat_P*Grid_para.A_b/1000; 0];
%PERUN
setpoint_PERUN = [solution.PV_P_perun*Grid_para.A_b/1000; 0];
%SOLARMAX
setpoint_SOLARMAX = [solution.PV_P_solarmax*Grid_para.A_b/1000; 0];
%TOTAL
SM_info.setpoints = [setpoint_AFE;setpoint_SUPERCAP;setpoint_BESS;setpoint_PERUN;setpoint_SOLARMAX]';
sendToStateMachine(SM_info,"ISLAND")

solution.Iabs(21)*Grid_para.I_b

solution.variables = variables;
solution.duration = toc(T_outer);
solution.mode = SM_info.mode;
solution.SC_soc_actual = SM_info.SC_soc;
solution.time = time_loop;


disp(['Control finished in ', num2str(round(solution.duration,2)), ' seconds'])

%% Save all data
name = ['data_',num2str(n_outer)];
assignin('base',name, solution)
save(FileName,name,'-append')
clear(genvarname(name))



time_to_wait = Ts - toc(T_outer);
pause(time_to_wait);



clear yalmip


end

