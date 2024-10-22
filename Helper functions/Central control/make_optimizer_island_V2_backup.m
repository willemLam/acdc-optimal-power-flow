function [solution,OPF_termination] = make_optimizer_island_V2(step,Grid_para,idx,constraints,variables,prev_solution,OPF_termination)

%% initialize yalmip
optimal.IC1ac_P = sdpvar(1,1); %B2
optimal.IC1dc_P = sdpvar(1,1); %B3
optimal.IC1ac_Q = sdpvar(1,1); %B2

optimal.IC2ac_P = sdpvar(1,1); %B2
optimal.IC2dc_P = sdpvar(1,1); %B3
optimal.IC2ac_Q = sdpvar(1,1); %B2

optimal.IC3ac_P = sdpvar(1,1); %B2
optimal.IC3dc_P = sdpvar(1,1); %B3
optimal.IC3ac_Q = sdpvar(1,1); %B2

optimal.IC4ac_P = sdpvar(1,1); %B2
optimal.IC4dc_P = sdpvar(1,1); %B3
optimal.IC4ac_Q = sdpvar(1,1); %B2

optimal.SC_P = sdpvar(1,1); %B4
optimal.Samsung_P = sdpvar(1,1); %B4
optimal.Perun_P = sdpvar(1,1); %B1

optimal.Eabs = sdpvar(Grid_para.n_nodes,1);
optimal.Eang = sdpvar(Grid_para.n_nodes,1);

optimal.Ploss = sdpvar(1,1);
optimal.Qloss = sdpvar(1,1);

%% Compute sensitivity coef using E
idxCtrl = 1:Grid_para.n_nodes;
[SC, ~] = SC_voltage_rectangular_V4(variables.E,idx,Grid_para,idxCtrl);

[K.Eabs.P, K.Eabs.Q, K.Eabs.Eabs,  K.Eabs.Eang, K.Eang.P, K.Eang.Q, K.Eang.Eabs, K.Eang.Eang] = transform_K_polar(SC, Grid_para, idx);
[K.E.P, K.E.Q, K.E.Eabs, K.E.Eang] = transform_K_complex(SC, Grid_para, idx);
[K.Iabs.P, K.I.P , K.Iabs.Q, K.I.Q , K.Iabs.Eabs, K.I.Eabs , K.Iabs.Eang, K.I.Eang  ] = Coeffs_Currents_GPT(variables.I,variables.E,K.E,Grid_para) ;
[K.Rloss.P, K.Rloss.Q, K.Rloss.Eabs, K.Rloss.Eang, K.Xloss.P, K.Xloss.Q,  K.Xloss.Eabs, K.Xloss.Eang ] = Coeffs_Losses(variables.E, K.E, Grid_para);

%% Constraints
con = [];

% Interfacing Converters
con = [con (optimal.IC1ac_P <= constraints.IC1_Qmax)&(optimal.IC1ac_P >= constraints.IC1_Qmin)]; %B19
con = [con (optimal.IC1ac_Q <= constraints.IC1_Qmax)&(optimal.IC1ac_Q >= constraints.IC1_Qmin)]; %B19
con = [con (optimal.IC1dc_P <= constraints.IC1_Qmax)&(optimal.IC1dc_P >= constraints.IC1_Qmin)]; %B19

con = [con (optimal.IC2ac_P <= constraints.IC2_Qmax)&(optimal.IC2ac_P >= constraints.IC2_Qmin)]; %B19
con = [con (optimal.IC2ac_Q <= constraints.IC2_Qmax)&(optimal.IC2ac_Q >= constraints.IC2_Qmin)]; %B19
con = [con (optimal.IC2dc_P <= constraints.IC2_Qmax)&(optimal.IC2dc_P >= constraints.IC2_Qmin)]; %B19

con = [con (optimal.IC3ac_P <= constraints.IC3_Qmax)&(optimal.IC3ac_P >= constraints.IC3_Qmin)]; %B19
con = [con (optimal.IC3ac_Q <= constraints.IC3_Qmax)&(optimal.IC3ac_Q >= constraints.IC3_Qmin)]; %B19
con = [con (optimal.IC3dc_P <= constraints.IC3_Qmax)&(optimal.IC3dc_P >= constraints.IC3_Qmin)]; %B19

con = [con (optimal.IC4ac_P <= constraints.IC4_Qmax)&(optimal.IC4ac_P >= constraints.IC4_Qmin)]; %B19
con = [con (optimal.IC4ac_Q <= constraints.IC4_Qmax)&(optimal.IC4ac_Q >= constraints.IC4_Qmin)]; %B19
con = [con (optimal.IC4dc_P <= constraints.IC4_Qmax)&(optimal.IC4dc_P >= constraints.IC4_Qmin)]; %B19


% Storage
con = [con (optimal.SC_P <= constraints.SupCap_Pmax)&(optimal.SC_P >= constraints.SupCap_Pmin)]; %B9
con = [con (optimal.Samsung_P <= constraints.Sambat_Pmax)&(optimal.Samsung_P >= constraints.Sambat_Pmin)]; %B9
con = [con (optimal.Perun_P <= variables.P(9))&(optimal.Perun_P >= 0)]; %B19


optimal.P = [zeros(1,1); 
             zeros(1,1);
             variables.P(3); 
             zeros(1,1);
             variables.P(5);
             zeros(1,1);
             optimal.Samsung_P;
             zeros(1,1);
             optimal.Perun_P;
             zeros(1,1);
             variables.P(11);
             zeros(1,1);
             variables.P(13);
             variables.P(14);
             zeros(4,1);
             zeros(4,1);
             zeros(4,1);
%              optimal.IC1ac_P;
%              optimal.IC2ac_P;
%              optimal.IC3ac_P;
%              optimal.IC4ac_P;
%              optimal.IC1dc_P;
%              optimal.IC2dc_P;
%              optimal.IC3dc_P;
%              optimal.IC4dc_P;
             optimal.SC_P;
             variables.P(28);
             variables.P(29);
             variables.P(30) 
             ] ;
                                            % DC nodes
optimal.Q = [zeros(1,1); 
             zeros(1,1);
             variables.Q(3); 
             zeros(1,1);
             variables.Q(5);
             zeros(1,1);
             zeros(1,1);
             zeros(1,1);
             variables.Q(9);
             zeros(1,1);
             variables.Q(11);
             zeros(1,1);
             variables.Q(13);
             variables.Q(14);
             zeros(4,1);
             optimal.IC1ac_Q;
             optimal.IC2ac_Q;
             optimal.IC3ac_Q;
             optimal.IC4ac_Q;
             zeros(8,1);
             ] ;

% con = [con optimal.SC_P == variables.Eabs(26) * (Grid_para.YY(22,26)*optimal.Eabs(22) + Grid_para.YY(26,26)*variables.Eabs(26) + Grid_para.YY(27,26)*optimal.Eabs(27)) ];

% con = [con optimal.IC1ac_P==-optimal.IC1dc_P];
% con = [con optimal.IC2ac_P==-optimal.IC2dc_P];
% con = [con optimal.IC3ac_P==-optimal.IC3dc_P];
% con = [con optimal.IC4ac_P==-optimal.IC4dc_P];

% Voltage magnitude constraints


        

con = [con optimal.Eabs == variables.Eabs - K.Eabs.P*(variables.P - optimal.P)  ...
                                          - K.Eabs.Q*(variables.Q - optimal.Q)  ...
                                          - K.Eabs.Eabs*(variables.Eabs  - optimal.Eabs)  ...
                                          - K.Eabs.Eang*(variables.Eang  - optimal.Eang)];
% K.Eabs.P(18,:)

con = [con optimal.Eabs(idx.vscac_vv) ==  1];
con = [con optimal.Eabs(idx.vdc) ==  0.9375];
% con = [con optimal.Eabs(idx.vscdc_vq) ==  0.9372 + 0*variables.Eabs(idx.vscdc_vq)];
             
con = [con ( optimal.Eabs <= constraints.E_max)];
con = [con ( optimal.Eabs >= constraints.E_min)];



% Voltage angle constraints
             
con = [con optimal.Eang == variables.Eang - K.Eang.P*(variables.P - optimal.P)  ...
                                          - K.Eang.Q*(variables.Q - optimal.Q)  ...
                                          - K.Eang.Eabs*(variables.Eabs  - optimal.Eabs)  ...
                                          - K.Eang.Eang*(variables.Eang  - optimal.Eang)]; %update non controllable V nodes

con = [con optimal.Eang(idx.vscac_vv) ==  0];
con = [con optimal.Eang(idx.vscdc_vv) ==  0];
con = [con optimal.Eang(idx.vdc) ==  0];
con = [con optimal.Eang(idx.pdc) ==  0];
con = [con optimal.Eang(idx.vscdc_vq) ==  0];            
             
con = [con (optimal.Eang <= pi/2)&(optimal.Eang >= -pi/2)]; 

% Current flow constraints
optimal.Iabs = variables.Iabs - K.Iabs.P*(variables.P - optimal.P)  ...
                              - K.Iabs.Q*(variables.Q - optimal.Q)  ...
                              - K.Iabs.Eabs*(variables.Eabs  - optimal.Eabs)  ...
                              - K.Iabs.Eang*(variables.Eang  - optimal.Eang); %update non controllable V nodes

con = [con (optimal.Iabs <= constraints.I_max)&(optimal.Iabs >= -0*constraints.I_max)];

% Losses for slack
variables.Ploss = real(sum(variables.S));
variables.Qloss = imag(sum(variables.S));

con = [con (optimal.Ploss == variables.Ploss -  K.Rloss.P*(variables.P - optimal.P)...
                                             -  K.Rloss.Q*(variables.Q - optimal.Q)...
                                             -  K.Rloss.Eabs*(variables.Eabs  - optimal.Eabs)...
                                             -  K.Rloss.Eang*(variables.Eang  - optimal.Eang))];

con = [con (optimal.Qloss == variables.Qloss -  K.Xloss.P*(variables.P - optimal.P)...
                                             -  K.Xloss.Q*(variables.Q - optimal.Q)...
                                             -  K.Xloss.Eabs*(variables.Eabs  - optimal.Eabs)...
                                             -  K.Xloss.Eang*(variables.Eang  - optimal.Eang))];

% con = [con optimal.Ploss <= 0];     %just to be sure

% con = [con optimal.SC_P == optimal.Eabs(27) * (Grid_para.G(23,27)*optimal.Eabs(23) + Grid_para.G(27,27)*optimal.Eabs(27) + Grid_para.G(28,27)*optimal.Eabs(28))];


con = [con optimal.SC_P == -(sum(optimal.P([1:26,28:30])) - optimal.Ploss)];
con = [con optimal.IC1ac_Q == -(sum(optimal.Q([1:18,20:22])) - optimal.Qloss)];

con = [con abs(variables.Eabs - optimal.Eabs) <= step];
con = [con abs(variables.Eang - optimal.Eang) <= step];
con = [con abs(variables.Iabs - optimal.Iabs) <= step];

% Extra constrain to avoid loop flows
con = [con sum((optimal.Iabs(18:21)))  >= 1.2*max(optimal.Iabs(18:21))] ;    % 1.2 as little margin for e.g. losses
con = [con sum((optimal.Iabs(22:25)))  >= 1.2*max(optimal.Iabs(22:25))] ;   
% con = [con sum((optimal.Iabs(18:21)))  >= max(optimal.Iabs(18:21))];

%% cost function

%minimize losses
cost.loss = norm(optimal.Ploss ,2) + norm(optimal.Qloss ,2); 

cost.IC_V = norm([0.9375 - optimal.Eabs(idx.vscdc_vq)]);
cost.IC_P = norm([optimal.Iabs(18:25)]);
cost.SC_P = norm([optimal.SC_P]);

cost.ICac_Q = norm([optimal.IC1ac_Q]);

cost.sum = cost.loss  + cost.ICac_Q +  cost.SC_P + cost.IC_P + 10*cost.IC_V;

% Solve the problem
options = sdpsettings('solver','gurobi','verbose', 2);
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

% optimal.Iabs = zeros(2,1);
% optimal.Eang  = zeros(4,1);
solution.P = double(optimal.P);
solution.Q = double(optimal.Q);
solution.S = complex(solution.P,solution.Q);
solution.Eabs = double(optimal.Eabs);
solution.Eang = double(optimal.Eang);
solution.E = solution.Eabs .* exp(1i*solution.Eang);
solution.Iabs = double(optimal.Iabs);

solution.Ploss = double(optimal.Ploss);
solution.Qloss = double(optimal.Qloss);
solution.P_sambat = double(optimal.SC_P);



tol = 1e-7; n_max = 20;
[Grid_para, Filter_para, idx1, idx3, constraints] = initialize_RT_island();
[solution.reconstructed.E,~,~] = NR_rectangularACDC_1ph_general_V2(Grid_para,Filter_para,solution.S,solution.E,solution.E,idx1,tol,n_max);
solution.reconstructed.S = solution.reconstructed.E.*conj(Grid_para.YY*solution.reconstructed.E);
solution.reconstructed.Iabs = abs(get_Current_flow(solution.reconstructed.E,Grid_para));





solution.reconstructed.E;
% solution.reconstructed.S
% solution.reconstructed.Iabs;
end
