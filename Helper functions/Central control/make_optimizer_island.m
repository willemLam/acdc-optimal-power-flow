function [solution,OPF_termination] = make_optimizer_island(Grid_para,idx,constraints,variables,prev_solution,OPF_termination)

%% initialize yalmip
optimal.slack.P = sdpvar(1,1); %B27
optimal.slack.Q = sdpvar(1,1); %B19

optimal.PV_P_perun = sdpvar(1,1); %B11 Perun
optimal.PV_Q_perun = sdpvar(1,1); %B11 Perun
optimal.PV_P_solarmax = sdpvar(1,1); %B11 SolarMax
optimal.PV_P_emul4 = sdpvar(1,1); %B26 Emulator
optimal.PV_P_emul5 = sdpvar(1,1); %B25 Emulator
optimal.PV_P_emul6 = sdpvar(1,1); %B24 Emulator
optimal.PV_P_emul7 = sdpvar(1,1); %B23 Emulator

optimal.IC1ac_P = sdpvar(1,1); %B2
optimal.IC1dc_P = sdpvar(1,1); %B3
optimal.IC2ac_P = sdpvar(1,1); %B2
optimal.IC2dc_P = sdpvar(1,1); %B3
optimal.IC3ac_P = sdpvar(1,1); %B2
optimal.IC3dc_P = sdpvar(1,1); %B3
optimal.IC4ac_P = sdpvar(1,1); %B2
optimal.IC4dc_P = sdpvar(1,1); %B3

optimal.IC2_Q = sdpvar(1,1); %B17
optimal.IC3_Q = sdpvar(1,1); %B16
optimal.IC4_Q = sdpvar(1,1); %B18

optimal.Sambat_P = sdpvar(1,1); %B9

optimal.IC2_E = sdpvar(1,1); %B20
optimal.IC3_E = sdpvar(1,1); %B21
optimal.IC4_E = sdpvar(1,1); %B22

optimal.Eabs = sdpvar(Grid_para.n_nodes,1);
optimal.Eang = sdpvar(Grid_para.n_nodes,1);
% optimal.P = sdpvar(Grid_para.n_nodes,1);
% optimal.Q = sdpvar(Grid_para.n_nodes,1);


%% Compute sensitivity coef using E
idxCtrl = 1:Grid_para.n_nodes;
[SC, ~] = SC_voltage_rectangular_V4(variables.E,idx,Grid_para,idxCtrl);

[K.Eabs.P, K.Eabs.Q, K.Eabs.Eabs,  K.Eabs.Eang, K.Eang.P, K.Eang.Q, K.Eang.Eabs, K.Eang.Eang] = transform_K_polar(SC, Grid_para, idx);
[K.E.P, K.E.Q, K.E.Eabs, K.E.Eang] = transform_K_complex(SC, Grid_para, idx);
[K.Iabs.P, K.I.P , K.Iabs.Q, K.I.Q , K.Iabs.Eabs, K.I.Eabs , K.Iabs.Eang, K.I.Eang  ] = Coeffs_Currents_GPT(variables.I,variables.E,K.E,Grid_para) ;
[K.Rloss.P, K.Rloss.Q, K.Rloss.Eabs, K.Rloss.Eang, K.Xloss.P, K.Xloss.Q,  K.Xloss.Eabs, K.Xloss.Eang ] = Coeffs_Losses(variables.E, K.E, Grid_para);

%% Constraints
con = [];

% Photo Voltaic 
con = [con (optimal.PV_P_perun <= variables.PVmax_P_perun)&(optimal.PV_P_perun >= 0)]; %B7
con = [con (optimal.PV_P_solarmax <= variables.PVmax_P_solarmax)&(optimal.PV_P_solarmax >= 0)]; %B11
con = [con (optimal.PV_P_emul4 <= variables.PVmax_P_emul4)&(optimal.PV_P_emul4 >= 0)]; %B26
con = [con (optimal.PV_P_emul5 <= variables.PVmax_P_emul5)&(optimal.PV_P_emul5 >= 0)]; %B25
con = [con (optimal.PV_P_emul6 <= variables.PVmax_P_emul6)&(optimal.PV_P_emul6 >= 0)]; %B24
con = [con (optimal.PV_P_emul7 <= variables.PVmax_P_emul7)&(optimal.PV_P_emul7 >= 0)]; %B23

% con = [con (optimal.PV_Q_perun <= constraints.PV_Q_perun_max)&(optimal.PV_Q_perun >= constraints.PV_Q_perun_min)]; %B7
con = [con (optimal.PV_Q_perun <= 0*optimal.PV_P_perun)&(optimal.PV_Q_perun >= -0*optimal.PV_P_perun)]; %B7

% Interfacing Converters
con = [con (optimal.slack.Q <= constraints.IC1_Qmax)&(optimal.slack.Q >= constraints.IC1_Qmin)]; %B19
con = [con (optimal.IC2_Q <= constraints.IC2_Qmax)&(optimal.IC2_Q >= constraints.IC2_Qmin)]; %B20
con = [con (optimal.IC3_Q <= constraints.IC3_Qmax)&(optimal.IC3_Q >= constraints.IC3_Qmin)]; %B21
con = [con (optimal.IC4_Q <= constraints.IC4_Qmax)&(optimal.IC4_Q >= constraints.IC4_Qmin)]; %B22

% Storage
con = [con (optimal.Sambat_P <= constraints.Sambat_Pmax)&(optimal.Sambat_P >= constraints.Sambat_Pmin)]; %B9
con = [con (optimal.slack.P  <= constraints.SupCap_Pmax)&(optimal.slack.P  >= constraints.SupCap_Pmin)]; %B27



% Nodal active power injections
% optimal.P(9) = optimal.PV_P_perun + optimal.Sambat_P;
% optimal.P(11) = optimal.PV_P_solarmax + variables.PVmax_P_facade; %facade is not curtailable
% optimal.P(30) = optimal.PV_P_emul4;
% optimal.P(29) = optimal.PV_P_emul5;
% optimal.P(28) = optimal.PV_P_emul6;
% optimal.P(27) = optimal.PV_P_emul7 + optimal.slack.P ;
% 
% optimal.P(1) = 0;                                               % GCP is a zero injection node
% optimal.P(3) = variables.P(3);                                  % uncontrolable load
% optimal.P(5) = variables.P(5);                                  % uncontrolable leclanche
% optimal.P(13) = variables.P(13);                                % uncontrolable FC
% optimal.P(14) = variables.P(14);                                % uncontrolable EV
% optimal.P([2,4,6,7,8,10,12]) = 0;                               % zero injection node
% optimal.P([15,16,17,18]) = variables.P([15,16,17,18]);          % filter nodes

optimal.P = [zeros(1,1); 
             zeros(1,1);
             variables.P(3);
             zeros(1,1);
             variables.P(5);
             zeros(3,1);
             optimal.PV_P_perun + optimal.Sambat_P; 
             zeros(1,1);
             optimal.PV_P_solarmax + variables.PVmax_P_facade; 
             zeros(1,1);
             variables.P(13); 
             variables.P(14);
             zeros(4,1);
             optimal.IC1ac_P;
             optimal.IC2ac_P;
             optimal.IC3ac_P;
             optimal.IC4ac_P;
             optimal.IC1dc_P;
             optimal.IC2dc_P;
             optimal.IC3dc_P;
             optimal.IC4dc_P;
             optimal.PV_P_emul7 + optimal.slack.P;
             optimal.PV_P_emul6;
             optimal.PV_P_emul5;
             optimal.PV_P_emul4
             ] ;
         

con = [con optimal.IC1ac_P==-optimal.IC1dc_P];
con = [con optimal.IC2ac_P==-optimal.IC2dc_P];
con = [con optimal.IC3ac_P==-optimal.IC3dc_P];
con = [con optimal.IC4ac_P==-optimal.IC4dc_P];

% Nodal reactive power injections
% optimal.Q(9) = optimal.PV_Q_perun;                              % Perun reactive
% optimal.Q(19) = optimal.slack.Q;                                % Q slack
% optimal.Q([20,21,22]) = [optimal.IC1_Q,optimal.IC2_Q ,optimal.IC3_Q]; % AFE in Vdc-Q mode
% 
% optimal.Q(1) = 0;                                               % GCP is a zero injection node
% optimal.Q(3) = variables.Q(3);                                  % uncontrolable load
% optimal.Q(5) = variables.Q(5);                                  % uncontrolable leclanche
% optimal.Q(11) = variables.Q(11);                                % uncontrolable PV
% optimal.Q(13) = variables.Q(13);                                % uncontrolable FC
% optimal.Q(14) = variables.Q(14);                                % uncontrolable EV
% optimal.Q([2,4,6,7,8,10,12]) = 0;                               % zero injection node
% optimal.Q([15,16,17,18]) = variables.Q([15,16,17,18]);          % filter nodes
% optimal.Q([23:30]) = 0;                                         % DC nodes

optimal.Q = [zeros(1,1); 
             zeros(1,1);
             variables.Q(3);
             zeros(1,1);
             variables.Q(5);
             zeros(3,1);
             optimal.PV_Q_perun; 
             zeros(1,1);
             variables.Q(11); 
             zeros(1,1);
             variables.Q(13); 
             variables.Q(14); 
             zeros(4,1);
             optimal.slack.Q;
             optimal.IC2_Q;
             optimal.IC3_Q;
             optimal.IC4_Q;
             zeros(8,1)
             ] ;


% Voltage magnitude constraints
optimal.Eabs(idx.vscdc_vq) = [optimal.IC2_E, optimal.IC3_E, optimal.IC4_E]; % AFE in Vdc-Q mode
optimal.Eabs(idx.vscac_vv) = variables.Eabs(idx.vscac_vv);%abs(variables.E_slack_sp);                     % AFE in V-V mode
optimal.Eabs(idx.vdc) = abs(variables.E_SupCap);                            % DC slack controling Vdc

optimal.Eabs = variables.Eabs - K.Eabs.P*(variables.P - optimal.P)  ...
                              - K.Eabs.Q*(variables.Q - optimal.Q)  ...
                              - K.Eabs.Eabs*(variables.Eabs  - optimal.Eabs)  ...
                              - K.Eabs.Eang*(variables.Eang  - optimal.Eang); %update non controllable V nodes

con = [con ( optimal.Eabs <= constraints.E_max)];
con = [con ( optimal.Eabs >= constraints.E_min)];


% Voltage angle constraints
% optimal.Eang([idx.vscdc_vq;idx.vscdc_pq;idx.pdc;idx.vdc]) = 0;              % DC nodes: zero angle
% optimal.Eang(idx.vscac_vv) = variables.Eang(idx.vscac_vv);%angle(variables.E_slack_sp);
% 
% optimal.Eang = variables.Eang - K.Eang.P*(variables.P - optimal.P)  ...
%                               - K.Eang.Q*(variables.Q - optimal.Q)  ...
%                               - K.Eang.Eabs*(variables.Eabs  - optimal.Eabs)  ...
%                               - K.Eang.Eang*(variables.Eang  - optimal.Eang); %update non controllable V nodes
% 
% con = [con (optimal.Eang <= pi/2)&(optimal.Eang >= -pi/2)]; 

% Current flow constraints
% optimal.Iabs = variables.Iabs - K.Iabs.P*(variables.P - optimal.P)  ...
%                               - K.Iabs.Q*(variables.Q - optimal.Q)  ...
%                               - K.Iabs.Eabs*(variables.Eabs  - optimal.Eabs)  ...
%                               - K.Iabs.Eang*(variables.Eang  - optimal.Eang); %update non controllable V nodes
% 
% con = [con (optimal.Iabs <= constraints.I_max)&(optimal.Iabs >= -constraints.I_max)];

% Losses for slack
variables.Ploss = real(sum(variables.S));
variables.Qloss = imag(sum(variables.S));

optimal.Ploss = variables.Ploss  -  K.Rloss.P*(variables.P - optimal.P)...
                                 -  K.Rloss.Q*(variables.Q - optimal.Q)...
                                 -  K.Rloss.Eabs*(variables.Eabs  - optimal.Eabs)...
                                 -  K.Rloss.Eang*(variables.Eang  - optimal.Eang);

optimal.Qloss = variables.Qloss  -  K.Xloss.P*(variables.P - optimal.P)...
                                 -  K.Xloss.Q*(variables.Q - optimal.Q)...
                                 -  K.Xloss.Eabs*(variables.Eabs  - optimal.Eabs)...
                                 -  K.Xloss.Eang*(variables.Eang  - optimal.Eang);

% con = [con optimal.Ploss <= 0];     %just to be sure

con = [con optimal.slack.P == -(sum(optimal.P(setdiff(1:Grid_para.n_nodes,idx.vdc))) + optimal.Ploss)];
con = [con optimal.slack.Q == -(sum(optimal.Q(setdiff(1:Grid_para.n_nodes,idx.vscac_vv))) + optimal.Qloss)];

% Extra constrain to avoid loop flows
% con = [con sum(abs(optimal.Iabs(22:25)))  <= max(optimal.Iabs(22:25))] ;    % 1.2 as little margin for e.g. losses
% con = [con sum(abs(optimal.Iabs(22:25)))  >= max(optimal.Iabs(22:25))];

%% cost function

%prevent E dc from going to high
Edc_nominal = variables.E_slack_sp; % 750/800; %
cost.IC.E = 10* (norm([ (optimal.IC2_E - Edc_nominal), (optimal.IC3_E - Edc_nominal), (optimal.IC4_E - Edc_nominal)],2));
% cost.Edc = norm(optimal.Eabs(23:30),2);

%min converters
% cost.IC.P =  norm(optimal.Iabs(22:25) ,2); %norm([optimal.P(idx.vscac_pq),optimal.P(idx.vscac_vv)],2);

%minimize losses
cost.loss = norm(optimal.Ploss ,2) + norm(optimal.Qloss ,2); 

%prevent curtailement
cost.PV_P =	norm([  (optimal.PV_P_perun -    variables.PVmax_P_perun     ); ...
                    (optimal.PV_P_solarmax - variables.PVmax_P_solarmax  ); ...
                    (optimal.PV_P_emul4 -    variables.PVmax_P_emul4     ); ...
                    (optimal.PV_P_emul5 -    variables.PVmax_P_emul5     ); ...
                    (optimal.PV_P_emul6 -    variables.PVmax_P_emul6     ); ...
                    (optimal.PV_P_emul7 -    variables.PVmax_P_emul7     )],2);

% not needed in theory (losses already accounts for them)
cost.IC.Q = norm([optimal.IC2_Q ,optimal.IC3_Q,optimal.IC4_Q],2);
cost.PV_Q =	norm(    optimal.PV_Q_perun ); 

cost.sambat = norm([optimal.Sambat_P]);

% minimize for island
cost.slack.P = norm( optimal.slack.P ); 
cost.slack.Q = norm( optimal.slack.Q ); 

% cost.sum = cost.IC.E + cost.loss  + cost.PV_P + cost.IC.Q + cost.PV_Q + cost.slack.Q; 
% cost.sum = cost.Edc + 0*cost.IC.E + cost.loss  + cost.slack.P + cost.slack.Q +  cost.PV_P + cost.IC.Q + cost.PV_Q + 0*cost.sambat; 
cost.sum = cost.IC.E + cost.loss +  cost.PV_P + cost.IC.Q; 
cost.sum = cost.IC.E + cost.loss  + cost.slack.Q +  cost.PV_P + cost.IC.Q + cost.sambat + cost.slack.P; 


% Solve the problem
options = sdpsettings('solver','gurobi','verbose', 0);
% options.DualReductions=0;
% options.InfUnbdInfo=1;


sol = optimize(con, cost.sum, options);


% fprintf('cost of total     %f \n', double(cost.sum))
% fprintf('cost of IC E    %f \n', double(cost.IC.E))
% fprintf('cost of losses     %f \n', double(cost.loss))
% fprintf('cost of curtailment     %f \n', double(cost.PV_P))
% fprintf('cost of IC Q     %f \n', double(cost.IC.Q))
% fprintf('cost of Samsung     %f \n', double(cost.sambat))
% fprintf('cost of slack P     %f \n', double(cost.slack.P))
% fprintf('cost of slack Q     %f \n', double(cost.slack.Q))

%%
if sol.problem == 0
 disp('Optimal solution found!');
 OPF_termination = 0;
elseif sol.problem ~= 0
 sol.info
 yalmiperror(sol.problem)
 warning('No optimal solution found - previous solution is implemented')
 solution = prev_solution;
 OPF_termination = OPF_termination+1;
end

optimal.Iabs = zeros(28,1);
optimal.Eang  = zeros(30,1);
solution.P = double(optimal.P);
solution.Q = double(optimal.Q);
solution.S = complex(solution.P,solution.Q);
solution.Eabs = double(optimal.Eabs);
solution.Eang = double(optimal.Eang);
solution.E = solution.Eabs .* exp(1i*solution.Eang);
solution.Iabs = double(optimal.Iabs);
solution.PVmax_P_perun = variables.PVmax_P_perun;
solution.PVmax_P_solarmax = variables.PVmax_P_solarmax;
solution.PVmax_P_facade = variables.PVmax_P_facade;
solution.PVmax_P_emul4 = variables.PVmax_P_emul4;
solution.PVmax_P_emul5 = variables.PVmax_P_emul5;
solution.PVmax_P_emul6 = variables.PVmax_P_emul6;
solution.PVmax_P_emul7 = variables.PVmax_P_emul7;
solution.PV_P_perun = double(optimal.PV_P_perun);
solution.PV_P_solarmax = double(optimal.PV_P_solarmax);
solution.PV_P_emul4 = double(optimal.PV_P_emul4);
solution.PV_P_emul5 = double(optimal.PV_P_emul5);
solution.PV_P_emul6 = double(optimal.PV_P_emul6);
solution.PV_P_emul7 = double(optimal.PV_P_emul7);
solution.Ploss = double(optimal.Ploss);
solution.Qloss = double(optimal.Qloss);
solution.P_sambat = double(optimal.Sambat_P);



tol = 1e-7; n_max = 20;
[Grid_para, Filter_para, idx1, idx3, constraints] = initialize_RT_island();
idx = idx1; %single phase equivalent
[solution.reconstructed.E,~,~] = NR_rectangularACDC_1ph_general_V2(Grid_para,Filter_para,solution.S,solution.E,solution.E,idx1,tol,n_max);
solution.reconstructed.S = solution.reconstructed.E.*conj(Grid_para.YY*solution.reconstructed.E);
solution.reconstructed.Iabs = abs(get_Current_flow(solution.reconstructed.E,Grid_para));

abs(solution.reconstructed.E)
end
