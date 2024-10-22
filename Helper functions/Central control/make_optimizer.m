function [solution,OPF_termination] = make_optimizer(Grid_para,idx,constraints,variables,prev_solution,OPF_termination)

% initialize yalmip
optimal.slack.P = sdpvar(1,1);
optimal.slack.Q = sdpvar(1,1);

optimal.PV_P_perun = sdpvar(1,1); %B11 Perun
optimal.PV_P_solarmax = sdpvar(1,1); %B11 SolarMax
optimal.PV_P_emul4 = sdpvar(1,1); %B26 Emulator
optimal.PV_P_emul5 = sdpvar(1,1); %B25 Emulator
optimal.PV_P_emul6 = sdpvar(1,1); %B24 Emulator
optimal.PV_P_emul7 = sdpvar(1,1); %B23 Emulat
optimal.PV_Q_perun = sdpvar(1,1); %B11 Perunor

optimal.IC1_Q = sdpvar(1,1); %B18
optimal.IC2_Q = sdpvar(1,1); %B17
optimal.IC3_Q = sdpvar(1,1); %B16
optimal.IC4_Q = sdpvar(1,1); %B15
% Q_opt_PV1_SolarMax = 0; %sdpvar(1,1, 'full'); %B11  SolarMax

optimal.IC1_E = sdpvar(1,1); %B22
optimal.IC2_E = sdpvar(1,1); %B21
optimal.IC3_E = sdpvar(1,1); %B20
optimal.IC4_E = sdpvar(1,1); %B19

optimal.Eabs = sdpvar(Grid_para.n_nodes,1);
optimal.Eang = sdpvar(Grid_para.n_nodes,1);


% Compute sensitivity coef using E
idxCtrl = 1:Grid_para.n_nodes;
[SC, ~] = SC_voltage_rectangular(variables.E,idx,Grid_para,idxCtrl);
[K.Eabs.P, K.Eabs.Q, K.Eabs.E,K.Eang.P, K.Eang.Q, K.Eang.E] = transform_K_polar(SC, Grid_para, idx);
[K.E.P, K.E.Q, K.E.E] = transform_K_complex(SC, Grid_para, idx);
[K.Iabs.P, ~ , K.Iabs.Q, ~ , K.Iabs.E, ~ ] = Coeffs_Currents_GPT(variables,K.E,Grid_para) ;
[C.Rloss.P , C.Rloss.Q, C.Rloss.E , C.Xloss.P , C.Xloss.Q, C.Xloss.E ] = Coeffs_Losses(variables.E, SC, Grid_para, idx);

% Constraints
con = [];
% Photo Voltaic 
%P
con = [con (optimal.PV_P_perun <= variables.PVmax_P_perun)&(optimal.PV_P_perun >= 0)]; %B7
con = [con (optimal.PV_P_solarmax <= variables.PVmax_P_solarmax)&(optimal.PV_P_solarmax >= 0)]; %B11



con = [con (optimal.PV_P_emul4 <= variables.PVmax_P_emul4)&(optimal.PV_P_emul4 >= 0)]; %B26
con = [con (optimal.PV_P_emul5 <= variables.PVmax_P_emul5)&(optimal.PV_P_emul5 >= 0)]; %B25
con = [con (optimal.PV_P_emul6 <= variables.PVmax_P_emul6)&(optimal.PV_P_emul6 >= 0)]; %B24
con = [con (optimal.PV_P_emul7 <= variables.PVmax_P_emul7)&(optimal.PV_P_emul7 >= 0)]; %B23
%Q
% con = [con (optimal.PV_Q_perun <= constraints.PV_Q_perun_max)&(optimal.PV_Q_perun >= constraints.PV_Q_perun_min)]; %B7
con = [con (optimal.PV_Q_perun <= 0.3*optimal.PV_P_perun)&(optimal.PV_Q_perun >= -0.3*optimal.PV_P_perun)]; %B7


% Interfacing Converters
con = [con (optimal.IC1_Q <= constraints.IC1_Qmax)&(optimal.IC1_Q >= constraints.IC1_Qmin)]; %B22
con = [con (optimal.IC2_Q <= constraints.IC2_Qmax)&(optimal.IC2_Q >= constraints.IC2_Qmin)]; %B21
con = [con (optimal.IC3_Q <= constraints.IC3_Qmax)&(optimal.IC3_Q >= constraints.IC3_Qmin)]; %B20
con = [con (optimal.IC4_Q <= constraints.IC4_Qmax)&(optimal.IC4_Q >= constraints.IC4_Qmin)]; %B19


% Make the nodal power injections
optimal.P = sdpvar(Grid_para.n_nodes,1);
optimal.P(idx.slack) = optimal.slack.P;
optimal.P(9) = optimal.PV_P_perun;
optimal.P(11) = optimal.PV_P_solarmax + variables.PVmax_P_facade;
optimal.P(30) = optimal.PV_P_emul4;
optimal.P(29) = optimal.PV_P_emul5;
optimal.P(28) = optimal.PV_P_emul6;
optimal.P(27) = optimal.PV_P_emul7;
optimal.P([2:8,10,12:22,23:26]) = variables.P([2:8,10,12:22,23:26]);

%optimal.Q = variable.Q;
optimal.Q = sdpvar(Grid_para.n_nodes,1);
optimal.Q(idx.slack) = optimal.slack.Q;
optimal.Q(idx.vscac_vq) = [optimal.IC1_Q,optimal.IC2_Q ,optimal.IC3_Q ,optimal.IC4_Q];
optimal.Q(9) = optimal.PV_Q_perun;
optimal.Q([2:8,10:18,23:30]) = variables.Q([2:8,10:18,23:30]);


% Voltage constraints
% nodal voltages magnitude
optimal.Eabs(idx.vscdc_vq) = [optimal.IC1_E, optimal.IC2_E, optimal.IC3_E, optimal.IC4_E];

optimal.Eabs = variables.Eabs - K.Eabs.P*(variables.P - optimal.P)  ...
                             - K.Eabs.Q*(variables.Q - optimal.Q)  ...
                             - K.Eabs.E*(variables.Eabs  - optimal.Eabs); %update non controllable V nodes

con = [con ( optimal.Eabs <= constraints.E_max)];
con = [con ( optimal.Eabs >= constraints.E_min)];


% nodal voltages angle
optimal.Eang([idx.vscdc_vq;idx.vscdc_pq;idx.pdc]) = 0;
optimal.Eang = variables.Eang - K.Eang.P*(variables.P - optimal.P)  ...
                             - K.Eang.Q*(variables.Q - optimal.Q)  ...
                             - K.Eang.E*(variables.Eabs  - optimal.Eabs); %update non controllable V nodes

con = [con (optimal.Eang <= pi/2)&(optimal.Eang >= -pi/2)]; 




% Slack constraints
con = [con (optimal.slack.P <= 1)&(optimal.slack.P >= -1)]; 
con = [con (optimal.slack.Q <= 1)&(optimal.slack.Q >= -1)]; 

% Current constraints
optimal.Iabs = variables.Iabs  - K.Iabs.P*(variables.P - optimal.P)  ...
                              - K.Iabs.Q*(variables.Q - optimal.Q)  ...
                              - K.Iabs.E*(variables.Eabs  - optimal.Eabs); %update non controllable V nodes

con = [con (optimal.Iabs <= constraints.I_max)&(optimal.Iabs >= -constraints.I_max)];

% losses

variables.Ploss = real(sum(variables.S));
variables.Qloss = imag(sum(variables.S));


optimal.Ploss = variables.Ploss  +  C.Rloss.P*(variables.P - optimal.P)...
                                +  C.Rloss.Q*(variables.Q - optimal.Q)...
                                +  C.Rloss.E*(variables.Eabs  - optimal.Eabs);

optimal.Qloss = variables.Qloss  +  C.Xloss.P*(variables.P - optimal.P)...
                                +  C.Xloss.Q*(variables.Q - optimal.Q)...
                                +  C.Xloss.E*(variables.Eabs  - optimal.Eabs);

% con = [con optimal.Ploss <= 0];


% slack
con = [con optimal.slack.P == -(sum(optimal.P(2:end)) - optimal.Ploss)];
con = [con optimal.slack.Q == -(sum(optimal.Q(2:end)) - optimal.Qloss)];


% cost function

%prevent E dc from going to high
Edc_nominal = 750/800;
cost.IC.E = norm([(optimal.IC1_E - Edc_nominal), (optimal.IC2_E - Edc_nominal), (optimal.IC3_E - Edc_nominal), (optimal.IC4_E - Edc_nominal)],2);

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
cost.IC.Q = norm([optimal.IC1_Q,optimal.IC2_Q ,optimal.IC3_Q ,optimal.IC4_Q],2);
cost.PV_Q =	norm(    optimal.PV_Q_perun ); 
cost.slack.Q = norm( optimal.slack.Q ); 

% cost.sum = cost.IC.E + cost.loss  + cost.PV_P + cost.IC.Q + cost.PV_Q + cost.slack.Q; 
cost.sum = cost.IC.E + 0*cost.loss  + cost.PV_P + cost.IC.Q + cost.PV_Q + cost.slack.Q; 


% Solve the problem
options = sdpsettings('solver','gurobi','verbose', 0);
% options.DualReductions=0;
% options.InfUnbdInfo=1;


sol = optimize(con, cost.sum, options);

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

end
