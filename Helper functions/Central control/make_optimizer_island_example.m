function [solution,OPF_termination] = make_optimizer_island_example(Grid_para,idx,constraints,variables,prev_solution,OPF_termination)

%% initialize yalmip

optimal.load_Q = sdpvar(1,1); %B1
optimal.PV1 = sdpvar(1,1); %B1

optimal.IC1ac_P = sdpvar(1,1); %B2
optimal.IC1dc_P = sdpvar(1,1); %B3
optimal.IC1ac_Q = sdpvar(1,1); %B2

optimal.IC2ac_P = sdpvar(1,1); %B2
optimal.IC2dc_P = sdpvar(1,1); %B3
optimal.IC2ac_Q = sdpvar(1,1); %B2

optimal.SC_P = sdpvar(1,1); %B4

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
% Storage
con = [con (optimal.SC_P <= constraints.SupCap_Pmax)&(optimal.SC_P >= constraints.SupCap_Pmin)]; %B9

% PV
con = [con (optimal.PV1 <= variables.P(1))&(optimal.PV1 >= 0)]; %B19


optimal.P = [optimal.PV1; 
             zeros(1,1);
             variables.P(3);
             zeros(1,1);
             zeros(1,1);
             optimal.SC_P;
             zeros(1,1);
             variables.P(8) 
             ] ;
                                            % DC nodes
optimal.Q = [optimal.load_Q; 
             zeros(1,1);
             variables.Q(3);
             optimal.IC2ac_Q;
             zeros(1,1);
             zeros(1,1);
             zeros(1,1);
             zeros(1,1)
             ] ;

     
% Voltage magnitude constraints


        

con = [con optimal.Eabs == variables.Eabs - K.Eabs.P*(variables.P - optimal.P)  ...
                              - K.Eabs.Q*(variables.Q - optimal.Q)  ...
                              - K.Eabs.Eabs*(variables.Eabs  - optimal.Eabs)  ...
                              - K.Eabs.Eang*(variables.Eang  - optimal.Eang)];

con = [con optimal.Eabs(idx.vscac_vv) ==  variables.Eabs(idx.vscac_vv)];
con = [con optimal.Eabs(idx.vdc) ==  0.9372 + 0*variables.Eabs(idx.vdc)];
% con = [con optimal.Eabs(idx.vscdc_vq) ==  0.9372 + 0*variables.Eabs(idx.vscdc_vq)];
             
con = [con ( optimal.Eabs <= constraints.E_max)];
con = [con ( optimal.Eabs >= constraints.E_min)];

   
% con = [con optimal.IC1ac_P == optimal.Eabs(1) * variables.Eabs(2) * ( real(Grid_para.YY(1,2)) +  imag(Grid_para.YY(1,2)) * (variables.Eang(2) - optimal.Eang(1))) + variables.Eabs(2) * variables.Eabs(2) * real(Grid_para.YY(2,2))];
% con = [con optimal.IC1ac_Q == optimal.Eabs(1) * variables.Eabs(2) * (-imag(Grid_para.YY(1,2)) +  real(Grid_para.YY(1,2)) * (variables.Eang(2) - optimal.Eang(1))) - variables.Eabs(2) * variables.Eabs(2) * imag(Grid_para.YY(2,2))];

% con = [con optimal.SC_P == optimal.Eabs(6) * (Grid_para.G(5,6)*optimal.Eabs(5) + Grid_para.G(6,6)*optimal.Eabs(6) + Grid_para.G(8,6)*optimal.Eabs(8)) ];

% Voltage angle constraints
             
con = [con optimal.Eang == variables.Eang - K.Eang.P*(variables.P - optimal.P)  ...
                              - K.Eang.Q*(variables.Q - optimal.Q)  ...
                              - K.Eang.Eabs*(variables.Eabs  - optimal.Eabs)  ...
                              - K.Eang.Eang*(variables.Eang  - optimal.Eang)]; %update non controllable V nodes

con = [con optimal.Eang(idx.vscac_vv) ==  variables.Eang(idx.vscac_vv)];
con = [con optimal.Eang(idx.vscdc_vv) ==  zeros(1,1)];
con = [con optimal.Eang(idx.vdc) ==  zeros(1,1)];
con = [con optimal.Eang(idx.vscdc_vq) ==  zeros(1,1)];            
             
con = [con (optimal.Eang <= pi/2)&(optimal.Eang >= -pi/2)]; 

% con = [con optimal.IC2ac_Q == variables.Q(4) ];
% con = [con optimal.IC1ac_Q == variables.Q(2) ];
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

con = [con (optimal.Ploss == variables.Ploss  -  K.Rloss.P*(variables.P - optimal.P)...
                                 -  K.Rloss.Q*(variables.Q - optimal.Q)...
                                 -  K.Rloss.Eabs*(variables.Eabs  - optimal.Eabs)...
                                 -  K.Rloss.Eang*(variables.Eang  - optimal.Eang))];

con = [con (optimal.Qloss == variables.Qloss  -  K.Xloss.P*(variables.P - optimal.P)...
                                 -  K.Xloss.Q*(variables.Q - optimal.Q)...
                                 -  K.Xloss.Eabs*(variables.Eabs  - optimal.Eabs)...
                                 -  K.Xloss.Eang*(variables.Eang  - optimal.Eang))];

% con = [con optimal.Ploss <= 0];     %just to be sure

con = [con optimal.SC_P == -(sum(optimal.P([1,2,3,4,5,7,8])) - optimal.Ploss)];
con = [con optimal.IC1ac_Q == -(sum(optimal.Q([1,3,4])) - optimal.Qloss)];

% Extra constrain to avoid loop flows
% con = [con sum(abs(optimal.Iabs(22:25)))  <= max(optimal.Iabs(22:25))] ;    % 1.2 as little margin for e.g. losses
% con = [con sum(abs(optimal.Iabs(22:25)))  >= max(optimal.Iabs(22:25))];

%% cost function

%minimize losses
cost.loss = norm(optimal.Ploss ,2) + norm(optimal.Qloss ,2); 

cost.SC_P = norm([optimal.SC_P]);

cost.ICac_Q = norm([optimal.IC1ac_Q]);

cost.sum = cost.loss  + cost.ICac_Q +  cost.SC_P ;

% Solve the problem
options = sdpsettings('solver','gurobi','verbose',0);
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

optimal.Iabs = zeros(2,1);
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
[Grid_para, Filter_para, idx1, idx3, constraints] = initialize_RT_island_example();
idx = idx1; %single phase equivalent
[solution.reconstructed.E,~,~] = NR_rectangularACDC_1ph_general_V2(Grid_para,Filter_para,solution.S,solution.E,solution.E,idx1,tol,n_max);
solution.reconstructed.S = solution.reconstructed.E.*conj(Grid_para.YY*solution.reconstructed.E);
solution.reconstructed.Iabs = abs(get_Current_flow(solution.reconstructed.E,Grid_para));
solution.reconstructed.E
solution.reconstructed.S
solution.reconstructed.Iabs

end
