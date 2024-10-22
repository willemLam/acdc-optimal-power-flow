function [solution,OPF_termination] = make_optimizer_connected_V3(Grid_para,Filter_para,idx,constraints,variables,prev_solution,OPF_termination,SM_info,Ts,state)

%% initialize yalmip
optimal.P_pcc = sdpvar(1,1); %B1
optimal.Q_pcc = sdpvar(1,1); %B1

% optimal.P_pv_perun = sdpvar(1,1); %B11 Perun
% optimal.PV_Q_perun = sdpvar(1,1); %B11 Perun
% optimal.P_pv_solarmax = sdpvar(1,1); %B11 SolarMax
% optimal.PV_P_emul4 = sdpvar(1,1); %B26 Emulator
% optimal.PV_P_emul5 = sdpvar(1,1); %B25 Emulator
% optimal.PV_P_emul6 = sdpvar(1,1); %B24 Emulator
% optimal.PV_P_emul7 = sdpvar(1,1); %B23 Emulator

optimal.P_quad_ic = sdpvar(3,1); %losses
optimal.Q_ic1 = sdpvar(1,1); %B18
optimal.Q_ic2 = sdpvar(1,1); %B17
optimal.Q_ic3 = sdpvar(1,1); %B16
%optimal.Q_ic4 = 0;%sdpvar(1,1); %B15

optimal.P_samsung = sdpvar(1,1); %B9
optimal.P_sc = sdpvar(1,1); %B27
% optimal.E_sc = sdpvar(1,1); %B27
optimal.SOC_sc = sdpvar(1,1); %B27

optimal.E_ic1 = sdpvar(1,1); %B22
% optimal.P_ic1 = 0; %B21
optimal.Pac_ic1 = sdpvar(1,1); %B21
optimal.Pac_ic2 = sdpvar(1,1); %B21
optimal.Pac_ic3 = sdpvar(1,1); %B20

optimal.Pdc_ic1 = sdpvar(1,1); %B21
optimal.Pdc_ic2 = sdpvar(1,1); %B21
optimal.Pdc_ic3 = sdpvar(1,1); %B20

% optimal.E_ic4 = sdpvar(1,1); %B19

optimal.Eabs = sdpvar(Grid_para.n_nodes,1);
optimal.Eang = sdpvar(Grid_para.n_nodes,1);
optimal.Iabs = sdpvar(Grid_para.n_lines,1);
optimal.Ploss = sdpvar(1,1); %B26
optimal.Qloss = sdpvar(1,1); %B26

%% Compute sensitivity coef using E
idxCtrl = 1:Grid_para.n_nodes;
[SC, ~] = SC_voltage_rectangular_V4(variables.E,idx,Grid_para,Filter_para,idxCtrl);


[K.Eabs.P, K.Eabs.Q, K.Eabs.Eabs,  K.Eabs.Eang, K.Eang.P, K.Eang.Q, K.Eang.Eabs, K.Eang.Eang] = transform_K_polar(SC, Grid_para, idx);
[K.E.P, K.E.Q, K.E.Eabs, K.E.Eang] = transform_K_complex(SC, Grid_para, idx);
[K.Iabs.P, K.I.P , K.Iabs.Q, K.I.Q , K.Iabs.Eabs, K.I.Eabs , K.Iabs.Eang, K.I.Eang  ] = Coeffs_Currents_GPT(variables.I,variables.E,K.E,Grid_para) ;
[K.Rloss.P, K.Rloss.Q, K.Rloss.Eabs, K.Rloss.Eang, K.Xloss.P, K.Xloss.Q,  K.Xloss.Eabs, K.Xloss.Eang ] = Coeffs_Losses(variables.E, K.E, Grid_para);

%% Constraints
con = [];

% Photo Voltaic 
% con = [con (optimal.P_pv_perun <= variables.PVmax_P_perun)&(optimal.P_pv_perun >= 0)]; %B7
% con = [con (optimal.P_pv_solarmax <= variables.PVmax_P_solarmax)&(optimal.P_pv_solarmax >= 0)]; %B11
% con = [con (optimal.P_pv_solarmax <= variables.PVmax_P_solarmax)&(optimal.P_pv_solarmax >= variables.PVmax_P_solarmax)]; %B11


% con = [con (optimal.PV_P_emul4 <= variables.PVmax_P_emul4)&(optimal.PV_P_emul4 >= 0)]; %B26
% con = [con (optimal.PV_P_emul5 <= variables.PVmax_P_emul5)&(optimal.PV_P_emul5 >= 0)]; %B25
% con = [con (optimal.PV_P_emul6 <= variables.PVmax_P_emul6)&(optimal.PV_P_emul6 >= 0)]; %B24
% con = [con (optimal.PV_P_emul7 <= variables.PVmax_P_emul7)&(optimal.PV_P_emul7 >= 0)]; %B23

% Interfacing Converters
con = [con (optimal.Q_ic1 <= constraints.IC1_Qmax)&(optimal.Q_ic1 >= constraints.IC1_Qmin)]; %B22
con = [con (optimal.Q_ic2 <= constraints.IC2_Qmax)&(optimal.Q_ic2 >= constraints.IC2_Qmin)]; %B21
con = [con (optimal.Q_ic3 <= constraints.IC3_Qmax)&(optimal.Q_ic3 >= constraints.IC3_Qmin)]; %B20
% con = [con (optimal.Q_ic4 <= constraints.IC4_Qmax)&(optimal.Q_ic4 >= constraints.IC4_Qmin)]; %B19

con = [con (optimal.Pac_ic1 <= constraints.IC1_Qmax)&(optimal.Pac_ic1 >= constraints.IC1_Qmin)]; %B20
con = [con (optimal.Pac_ic2 <= constraints.IC2_Qmax)&(optimal.Pac_ic2 >= constraints.IC2_Qmin)]; %B21
con = [con (optimal.Pac_ic3 <= constraints.IC3_Qmax)&(optimal.Pac_ic3 >= constraints.IC3_Qmin)]; %B20

con = [con (optimal.Pdc_ic1 <= constraints.IC1_Qmax)&(optimal.Pdc_ic1 >= constraints.IC1_Qmin)]; %B21
con = [con (optimal.Pdc_ic2 <= constraints.IC2_Qmax)&(optimal.Pdc_ic2 >= constraints.IC2_Qmin)]; %B21
con = [con (optimal.Pdc_ic3 <= constraints.IC3_Qmax)&(optimal.Pdc_ic3 >= constraints.IC3_Qmin)]; %B20


% Storage
% SC is 0.196875kwh
% 90% efficiency
con = [con (optimal.SOC_sc == SM_info.SC_soc/100 - 0.90*(optimal.P_sc*Ts/3600)/(196.875/Grid_para.A_b))]; %sign: + is inject, - is charge
con = [con (optimal.P_samsung <= constraints.Sambat_Pmax)&(optimal.P_samsung >= constraints.Sambat_Pmin)]; %B9

% supercap limits
Delta_SOC = abs(50 - SM_info.SC_soc);
% SupCap_Pmin_adaptive = max( min(-Delta_SOC * 50/Grid_para.A_b,(-Delta_SOC + 5) * 150/Grid_para.A_b),constraints.SupCap_Pmin);
% SupCap_Pmax_adaptive = min( max( Delta_SOC * 50/Grid_para.A_b,( Delta_SOC - 5) * 150/Grid_para.A_b),constraints.SupCap_Pmax);
SupCap_Pmin_adaptive = max(-Delta_SOC * 150/Grid_para.A_b,constraints.SupCap_Pmin);
SupCap_Pmax_adaptive = min(Delta_SOC * 150/Grid_para.A_b,constraints.SupCap_Pmax);


con = [con (optimal.P_sc <= SupCap_Pmax_adaptive)&(optimal.P_sc >= SupCap_Pmin_adaptive)]; %B27
con = [con (optimal.SOC_sc <= 0.9)&(optimal.SOC_sc >= 0.1)]; %B27



% Nodal power injection
optimal.P = [optimal.P_pcc; %PCC
             zeros(1,1);
             variables.P(3); %uncontrolable load
             zeros(1,1);
             variables.P(5); %uncontrolable battery
             zeros(3,1);
             -0*constraints.BESS_loss + optimal.P_samsung; %controlable PV and battery
             zeros(1,1);
             variables.P(11); %variables.PV_P_facade + optimal.P_pv_solarmax; %un/controlable PV
             zeros(1,1);
             variables.P(13); %uncontrolable FC and EL
             variables.P(14); %uncontrolable EVCS
             zeros(4,1);
             optimal.Pac_ic1; %non controlable
             optimal.Pac_ic2; %controllable IC2
             optimal.Pac_ic3; %controllable IC3
             optimal.Pdc_ic1; %non controlable
             optimal.Pdc_ic2; %controllable IC2
             optimal.Pdc_ic3; %controllable IC3
             optimal.P_sc; %controlable supercap
             variables.P(26);  %uncontrolable PV
             variables.P(27);  %uncontrolable PV
%              variables.P(30)   %uncontrolable PV
             ] ;



optimal.Q = [optimal.Q_pcc; %PCC
             zeros(1,1);
             variables.Q(3); %uncontrolable load
             zeros(1,1);
             variables.Q(5); %uncontrolable battery
             zeros(3,1);
             variables.Q(9); %uncontrolable PV and battery
             zeros(1,1);
             variables.Q(11); %uncontrolable PV
             zeros(1,1);
             variables.Q(13);  %uncontrolable FC and EL
             variables.Q(14);  %uncontrolable EVCS
             zeros(4,1);
             optimal.Q_ic1;  %controlable AFE1
             optimal.Q_ic2;  %controlable AFE2
             optimal.Q_ic3;  %controlable AFE3
%              optimal.Q_ic4;  %controlable AFE4
             zeros(6,1)
             ] ;



% Voltage magnitude constraints
con = [con (optimal.Eabs(idx.vscdc_vq) == [optimal.E_ic1])]; % AFE in Vdc-Q mode
con = [con (optimal.Eabs(idx.slack) == variables.Eabs(idx.slack))];

con = [con (optimal.Eabs == variables.Eabs - K.Eabs.P*(variables.P - optimal.P)  ...
                              - K.Eabs.Q*(variables.Q - optimal.Q)  ...
                              - K.Eabs.Eabs*(variables.Eabs  - optimal.Eabs)  ...
                              - 0*K.Eabs.Eang*(variables.Eang  - optimal.Eang))]; %update non controllable V nodes

con = [con ( optimal.Eabs <= constraints.E_max)];
con = [con ( optimal.Eabs >= constraints.E_min)];

% if SM_info.mode == "Prepare Island"
con = [con (optimal.Eabs(25) == constraints.Edc_nominal)];
% con = [con (optimal.Eabs(22) == constraints.Edc_nominal)];
% end

% Voltage angle constraints
con = [con (optimal.Eang == variables.Eang - K.Eang.P*(variables.P - optimal.P)  ...
                              - K.Eang.Q*(variables.Q - optimal.Q)  ...
                              - K.Eang.Eabs*(variables.Eabs  - optimal.Eabs)  ...
                              - 0*K.Eang.Eang*(variables.Eang  - optimal.Eang))]; %update non controllable V nodes

% Current flow constraints
con = [con (optimal.Iabs == variables.Iabs - K.Iabs.P*(variables.P - optimal.P)  ...
                              - K.Iabs.Q*(variables.Q - optimal.Q)  ...
                              - K.Iabs.Eabs*(variables.Eabs  - optimal.Eabs)  ...
                              - 0*K.Iabs.Eang*(variables.Eang  - optimal.Eang))];
                          
con = [con (optimal.Iabs <= constraints.I_max)&(optimal.Iabs >= -constraints.I_max)];

if Filter_para.Include_losses == 1
     Iac = abs(Grid_para.YY(Grid_para.pos_ac3(:,1),:) * variables.E);
     Idc = abs(Grid_para.YY(Grid_para.pos_dc3(:,1),:) * variables.E);
     Vdc = abs(variables.E(Grid_para.pos_dc3(:,1)));
    
    optimal.P_quad_ic = Filter_para.a .* Vdc + Filter_para.c .* Iac.^2 ...
                           + Filter_para.d .* Idc + Filter_para.e .* Idc.^2;
else
    optimal.P_quad_ic = zeros(3,1);
end

con = [con (optimal.Pdc_ic1 == -optimal.Pac_ic1 - optimal.P_quad_ic(1) )];
con = [con (optimal.Pdc_ic2 == -optimal.Pac_ic2 - optimal.P_quad_ic(2) )];
con = [con (optimal.Pdc_ic3 == -optimal.Pac_ic3 - optimal.P_quad_ic(3) )];


% Losses for slack
variables.Ploss = real(sum(variables.S));
variables.Qloss = imag(sum(variables.S));

con = [con (optimal.Ploss == variables.Ploss   - K.Rloss.P*(variables.P - optimal.P)...
                                  - K.Rloss.Q*(variables.Q - optimal.Q)...
                                  - K.Rloss.Eabs*(variables.Eabs  - optimal.Eabs)...
                                  - 0*K.Rloss.Eang*(variables.Eang  - optimal.Eang))];

con = [con (optimal.Qloss == variables.Qloss   - K.Xloss.P*(variables.P - optimal.P)...
                                  - K.Xloss.Q*(variables.Q - optimal.Q)...
                                  - K.Xloss.Eabs*(variables.Eabs  - optimal.Eabs)...
                                  - 0*K.Xloss.Eang*(variables.Eang  - optimal.Eang))];

% con = [con optimal.Ploss <= 0];     %just to be sure

con = [con optimal.P_pcc == (-sum(optimal.P(2:Grid_para.n_nodes)) + (optimal.Ploss))];
con = [con optimal.Q_pcc == (-sum(optimal.Q(2:Grid_para.n_nodes)) + (optimal.Qloss))];

% % Extra constrain to avoid loop flows
% con = [con sum(abs(optimal.Iabs(22:25)))  <= max(optimal.Iabs(22:25))] ;    % 1.2 as little margin for e.g. losses
% con = [con sum(abs(optimal.Iabs(22:25)))  >= -max(optimal.Iabs(22:25))];


%% extra constraint to limit rate only on line 22-25
% Delta_I = 1/Grid_para.I_b;
% con = [con (optimal.Iabs(21) <= state.Iabs(21)+Delta_I)&(optimal.Iabs(21) >= state.Iabs(21)-Delta_I)]; %B22

% Delta_Q = 500/Grid_para.A_b;
% con = [con (optimal.Q_ic1 <= imag(state.S(19))+Delta_Q)&(optimal.Q_ic1 >= imag(state.S(19))-Delta_Q)]; %B22
% con = [con (optimal.Q_ic2 <= imag(state.S(20))+Delta_Q)&(optimal.Q_ic1 >= imag(state.S(20))-Delta_Q)]; %B22
% con = [con (optimal.Q_ic3 <= imag(state.S(21))+Delta_Q)&(optimal.Q_ic1 >= imag(state.S(21))-Delta_Q)]; %B22


if SM_info.mode == "Prepare Island"
%      con = [con (optimal.Iabs(21) <= Delta_I)]; %B22
end

%% cost function

%prevent E dc from going to high
SOC_ref = 0.5;
% cost.E_afe = norm((optimal.E_ic1 - Edc_nominal),2)*0 + (norm([(optimal.E_ic2 - Edc_nominal), (optimal.E_ic3 - Edc_nominal), (optimal.E_ic4 - Edc_nominal)],2));

%% Voltages
E_min = 0;
E_max = 15/Grid_para.Vdc_b;

% cost.E_afe = norm((optimal.E_ic1 - constraints.Edc_nominal),2);
% cost.E_afe_norm = (cost.E_afe - E_min)/(E_max - E_min);

cost.E_dc = norm((optimal.Eabs(22:27) - constraints.Edc_nominal),2);
cost.E_dc_norm = (cost.E_dc - E_min)/(E_max - E_min);

% cost.E_sc = norm(optimal.Eabs(25) - constraints.Edc_nominal);
% cost.E_sc_norm =  (cost.E_sc - E_min)/(E_max - E_min);

% cost.I_afe = norm(2/Grid_para.I_b - optimal.Iabs(21),2);
cost.I_afe =  norm(optimal.Eabs(25) - optimal.Eabs(22));
cost.I_afe_norm = (cost.I_afe - E_min)/(E_max - E_min);

%% Supercap SOC
SOC_min = 0;
SOC_max = 0.2;

cost.SOC_sc = norm(SOC_ref - optimal.SOC_sc,2);
cost.SOC_sc_norm = (cost.SOC_sc - SOC_min)/(SOC_max - SOC_min);

%% losses
P_min = 0;
P_max = 10000/Grid_para.A_b;

cost.loss = norm([optimal.Ploss, optimal.Qloss] ,2); 
cost.loss_norm = (cost.loss - P_min)/(P_max - P_min);

%% Q AFE not needed in theory (losses already accounts for them)
cost.Q_afe = norm([optimal.Q_ic1,optimal.Q_ic2 ,optimal.Q_ic3],2);
cost.Q_afe23 = norm([optimal.Q_ic2 ,optimal.Q_ic3],2);
cost.Q_afe1 = norm([optimal.Q_ic1],2); %The current through the future afe should be 0
cost.Q_afe_norm =  (cost.Q_afe - P_min)/(P_max - P_min);
cost.Q_afe23_norm =  (cost.Q_afe23 - P_min)/(P_max - P_min);
cost.Q_afe1_norm =  (cost.Q_afe1 - P_min)/(P_max - P_min);

cost.P_afe23 = norm([optimal.Pac_ic2 ,optimal.Pac_ic3],2);
cost.P_afe23_norm = (cost.P_afe23 - P_min)/(P_max - P_min);

cost.P_samsung = norm(optimal.P_samsung,2);
cost.P_sc = norm(optimal.P_sc,2);
cost.P_samsung_norm = (cost.P_samsung - P_min)/(P_max - P_min);
cost.P_sc_norm = (cost.P_sc - P_min)/(P_max - P_min);

%% PCC
cost.PCC_p = norm( optimal.P_pcc,2 ); 
cost.PCC_q = norm( optimal.Q_pcc,2 ); 
cost.PCC_p_norm = (cost.PCC_p - P_min)/(P_max - P_min);
cost.PCC_q_norm =  (cost.PCC_q - P_min)/(P_max - P_min);

if     SM_info.mode == "Grid connected"
    cost.sum =   100*cost.I_afe_norm + 100*cost.SOC_sc_norm + 0.1*cost.E_dc_norm + 0*cost.loss_norm + 5*cost.Q_afe_norm + 0.5*cost.P_afe23_norm + 10*cost.PCC_q_norm + cost.P_samsung_norm; 
elseif SM_info.mode == "Prepare Island" % here
    cost.sum =   100*cost.I_afe_norm + 100*cost.SOC_sc_norm + 0.1*cost.E_dc_norm + 0*cost.loss_norm + 5*cost.Q_afe_norm + 0.5*cost.P_afe23_norm + 10*cost.PCC_q_norm + 0.1*cost.P_samsung_norm + cost.PCC_p_norm + 10*cost.Q_afe1_norm; 
elseif SM_info.mode == "Forming but connected" % here
    cost.sum =   100*cost.I_afe_norm + 100*cost.SOC_sc_norm + 0.1*cost.E_dc_norm + 0*cost.loss_norm + 5*cost.Q_afe_norm + 0.5*cost.P_afe23_norm + 10*cost.PCC_q_norm + 0.1*cost.P_samsung_norm + cost.PCC_p_norm + 10*cost.Q_afe1_norm; 
elseif SM_info.mode == "Island"
    cost.sum =   100*cost.I_afe_norm + 100*cost.SOC_sc_norm + 0.1*cost.E_dc_norm + 0*cost.loss_norm + 5*cost.Q_afe_norm + 0.5*cost.P_afe23_norm + 10*cost.PCC_q_norm + 0.1*cost.P_samsung_norm + cost.PCC_p_norm + 10*cost.Q_afe1_norm; 
elseif SM_info.mode == "Prepare for connection"
    cost.sum =   100*cost.I_afe_norm + 100*cost.SOC_sc_norm + 0.1*cost.E_dc_norm + 0*cost.loss_norm + 5*cost.Q_afe_norm + 0.5*cost.P_afe23_norm + 10*cost.PCC_q_norm + 0.1*cost.P_samsung_norm + cost.PCC_p_norm + 10*cost.Q_afe1_norm; 
else
    warning('Wrong operation mode')
end


% Solve the problem
options = sdpsettings('solver','gurobi','verbose', 0);
% options.DualReductions=0;
% options.InfUnbdInfo=1;


sol = optimize(con, cost.sum, options);

% fprintf('cost of total     %f \n', double(cost.sum))
% fprintf('cost of loss    %f \n', 0*double(cost.loss_norm))
% % fprintf('cost of E_afe    %f \n', 10*double(cost.E_afe))
% fprintf('cost of I_afe    %f \n', 0*double(cost.I_afe_norm))
% fprintf('cost of E_dc     %f \n', double(cost.E_dc_norm))
% fprintf('cost of SOC_sc     %f \n', 20*double(cost.SOC_sc_norm))
% fprintf('cost of  Q_afe23     %f \n', 5*double(cost. Q_afe_norm))
% % fprintf('cost of Q_afe2     %f \n', double(cost.Q_afe2))
% fprintf('cost of P_samsung     %f \n', double(cost.P_samsung_norm))
% fprintf('cost of PCC Q     %f \n', double(cost.PCC_q_norm))

%%
% %prevent E dc from going to high
% SOC_ref = 0.5;
% % cost.E_afe = norm((optimal.E_ic1 - Edc_nominal),2)*0 + (norm([(optimal.E_ic2 - Edc_nominal), (optimal.E_ic3 - Edc_nominal), (optimal.E_ic4 - Edc_nominal)],2));
% 
% cost.E_afe = norm((optimal.E_ic1 - constraints.Edc_nominal),2);
% 
% cost.E_sc = 20*norm(optimal.Eabs(25) - constraints.Edc_nominal);
% 
% 
% cost.I_afe = 10*norm(optimal.Eabs(25) - optimal.Eabs(22));
% % cost.I_afe = norm(2/Grid_para.I_b - optimal.Iabs(21),2);
% 
% %Supercap SOC
% cost.SOC_sc = norm(SOC_ref - optimal.SOC_sc,2);
% 
% %minimize losses
% cost.loss = norm([optimal.Ploss, optimal.Qloss] ,2); 
% 
% % not needed in theory (losses already accounts for them)
% cost.Q_afe = 10*norm([optimal.Q_ic1,optimal.Q_ic2 ,optimal.Q_ic3],2);
% cost. Q_afe23 = 10*norm([optimal.Q_ic2 ,optimal.Q_ic3],2);
% cost. Q_afe1 = 20*norm([optimal.Q_ic1],2); %The current through the future afe should be 0
% 
% cost.P_afe1 = norm([optimal.Pac_ic2 ,optimal.Pac_ic3],2);
% 
% cost.P_samsung = norm(optimal.P_samsung,2);
% cost.P_sc = norm(optimal.P_sc,2);
% 
% % minimize for island
% cost.PCC_p = 30*norm( optimal.P_pcc,2 ); 
% cost.PCC_q = 30*norm( optimal.Q_pcc,2 ); 
% 
% if     SM_info.mode == "Grid connected"
% %     cost.sum = 10*cost.SOC_sc + 10*cost.E_afe + cost.PCC_q; 
%     cost.sum =   1*cost.I_afe + 10*cost.SOC_sc + 10*cost.E_afe + 0.1*cost.loss + cost.P_afe1 + cost.Q_afe + cost.PCC_q + cost.P_samsung; 
% elseif SM_info.mode == "Prepare Island" % here
%     cost.sum =   1*cost.I_afe + 10*cost.SOC_sc + 10*cost.E_sc  + 0.1*cost.loss + cost.P_afe1 + cost. Q_afe23 + cost. Q_afe1 + cost.PCC_p + cost.PCC_q;
% elseif SM_info.mode == "Forming but connected" % here
%     cost.sum =   1*cost.I_afe + 10*cost.SOC_sc + 10*cost.E_afe + 0.1*cost.loss + cost.P_afe1 + cost. Q_afe23 + cost. Q_afe1 + cost.PCC_p + cost.PCC_q;
% elseif SM_info.mode == "Island"
%     cost.sum =   1*cost.I_afe + 10*cost.SOC_sc + 10*cost.E_afe + 0.1*cost.loss + cost.P_afe1 + cost. Q_afe23 + cost. Q_afe1 + cost.PCC_p + cost.PCC_q;
% elseif SM_info.mode == "Prepare for connection"
%     cost.sum =   1*cost.I_afe + 10*cost.SOC_sc + 10*cost.E_afe + 0.1*cost.loss + cost.P_afe1 + cost. Q_afe23 + cost. Q_afe1 + cost.PCC_q + cost.P_samsung; 
% else
%     warning('Wrong operation mode')
% end

%%
if sol.problem == 0
     disp('Optimal solution found!');
     OPF_termination = 0;

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
%     solution.PVmax_P_emul4 = variables.PVmax_P_emul4;
%     solution.PVmax_P_emul5 = variables.PVmax_P_emul5;
%     solution.PVmax_P_emul6 = variables.PVmax_P_emul6;
%     solution.PVmax_P_emul7 = variables.PVmax_P_emul7;
    solution.PV_P_perun = variables.PVmax_P_perun;
    solution.PV_P_solarmax = variables.PVmax_P_solarmax;
    solution.PV_P_facade = variables.PVmax_P_facade;
    solution.PV_P_emul4 = variables.P(25);
    solution.PV_P_emul5 = variables.P(26);
    solution.PV_P_emul6 = variables.P(27);
%     solution.PV_P_emul7 = variables.P(30);
    solution.Ploss = double(optimal.Ploss);
    solution.Qloss = double(optimal.Qloss);
    solution.Sambat_P = double(optimal.P_samsung);
    solution.SupCap_P = double(optimal.P_sc);
    solution.SupCap_SOC = double(optimal.SOC_sc);
    solution.IC_Q = double([optimal.Q_ic1; optimal.Q_ic2; optimal.Q_ic3]);
    solution.IC_E = double([optimal.E_ic1]);
    solution.IC_P = double([optimal.Pac_ic2; optimal.Pac_ic3]);
    solution.P_quad_ic = double(optimal.P_quad_ic);

elseif sol.problem ~= 0
     sol.info
     yalmiperror(sol.problem)
     fprintf(2,'No optimal solution found - previous solution is implemented \n')
     solution = prev_solution;
     OPF_termination = OPF_termination+1;
end



% tol = 1e-7; n_max = 100;
% [Grid_para, Filter_para, idx1, idx3, constraints] = initialize_RT();
% idx = idx1; 
% [solution.reconstructed.E,~,~] = NR_rectangularACDC_1ph_general_V2(Grid_para,Filter_para,solution.S,solution.E,solution.E,idx1,tol,n_max);
% solution.reconstructed.S = solution.reconstructed.E.*conj(Grid_para.YY*solution.reconstructed.E);
% solution.reconstructed.Iabs = abs(get_Current_flow(solution.reconstructed.E,Grid_para));
% 
% abs(solution.reconstructed.E)
end


