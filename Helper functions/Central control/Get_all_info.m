function variable = Get_all_info(TYPE,Grid_para,Filter_para,idx)

% Get PV


%         variable.PVmax_P_perun = 0;  %Perun is off
%         variable.PVmax_P_solarmax = 0; %previous PV (accounting for the curtailment)
%         variable.PVmax_P_facade = 0; %previous PV 
%         variable.PVmax_P_emul4 = 0;
%         variable.PVmax_P_emul5 = 0;
%         variable.PVmax_P_emul6 = 0;
%         variable.PVmax_P_emul7 = 0;

%% Get state of the grid

if TYPE == 'EXPERIMENT'
    port_SE = 35005;
    % [E, S, I] = Get_SE_data(TYPE,port_SE, Grid_para,i);
    [E, S, I] = get_states_from_SE(port_SE, Grid_para);
    E(abs(E)<0.5) = 0.9375;
       
    PVmax = Get_PV_data_RT(TYPE,Grid_para);
    variable.PVmax_P_perun = PVmax.P_perun;  %Perun is off
    variable.PVmax_P_solarmax = PVmax.P_solarmax; %previous PV (accounting for the curtailment)
    variable.PVmax_P_facade = PVmax.P_facade; %previous PV 
    variable.PVmax_P_emul4 = 0;
    variable.PVmax_P_emul5 = 0;
    variable.PVmax_P_emul6 = 0;
    variable.PVmax_P_emul7 = 0;

elseif TYPE == 'SIMULATION'
    [E,S,I] = get_fake_states(idx, Grid_para,Filter_para);
    
    variable.PVmax_P_perun = real(S(9));  %Perun is off
    variable.PVmax_P_solarmax = 0.7*real(S(11)); %previous PV (accounting for the curtailment)
    variable.PVmax_P_facade = 0.3*real(S(11)); %previous PV 
    variable.PVmax_P_emul4 = real(S(30));
    variable.PVmax_P_emul5 = real(S(29));
    variable.PVmax_P_emul6 = real(S(28));
    variable.PVmax_P_emul7 = 0;

end
variable.P = real(S);
variable.Q = imag(S);
variable.E = E;
variable.S = S;
variable.Eabs = abs(E);
variable.Eang = angle(E);
variable.Iabs = abs(I);
variable.I = I;


% tol = 1e-7; n_max = 100;
% [E_bus2,~,~] = NR_rectangularACDC_1ph_general(Grid_para,Filter_para,S_bus,E_bus,Grid_para.E_0,idx1,tol,n_max);
% S_bus2 = E_bus2.*conj(Grid_para.YY*E_bus2);
% I_line = get_Current_flow(E_bus,Grid_para);

%     Eabs = abs(E_bus);
%     Eang = angle(E_bus);
%     Eabs_slack = E_bus(idx1.slack);
% 
%     P_bus = real(S_bus);
%     Q_bus = imag(S_bus);
%     
%     I = abs(I_line);

%% Get data from Interfacing Converters
% [variable.IC_P,variable.IC_Q,variable.IC_Edc] = Get_IC_data(Grid_para);

end