function variable = Get_all_info_OFFLINE(TYPE,Grid_para,Filter_para,idx,n_outer,Ts,GHI,S_setpoint,E_setpoint)

% S_setpoint = prev_solution.S;
% E_setpoint = prev_solution.E;
% GHI = GHI_data(n_outer);

%% Get PV
    GI = GHI*1.1;
    GVI = GHI*1.2;
    Tair = 10;
    MatlabTime = datetime(2024,1,1,0,0,n_outer*Ts);
    UTC = 2;

    [PV_perun, PV_solarmax,PV_facade] = PV_production(GHI, GI, GVI, Tair, MatlabTime, UTC);

        variable.PVmax_P_perun = PV_perun/Grid_para.A_b; 
        variable.PVmax_P_solarmax = PV_solarmax/Grid_para.A_b; %previous PV (accounting for the curtailment)
        variable.PVmax_P_facade = PV_facade/Grid_para.A_b; %previous PV 
        variable.PVmax_P_emul4 = PV_perun/Grid_para.A_b/3;
        variable.PVmax_P_emul5 = PV_perun/Grid_para.A_b/3;
        variable.PVmax_P_emul6 = PV_perun/Grid_para.A_b/3;
        variable.PVmax_P_emul7 = PV_perun/Grid_para.A_b/3;

%% Get State
%% Validate solution
    tol = 1e-7; n_max = 100;
    [E,~,~] = NR_rectangularACDC_1ph_general(Grid_para,Filter_para,S_setpoint,E_setpoint,E_setpoint,idx,tol,n_max);
    S = E.*conj(Grid_para.YY*E);
    I = get_Current_flow(E,Grid_para);
        
        variable.P = real(S);
        variable.Q = imag(S);
        variable.E = E;
        variable.S = S;
        variable.Eabs = abs(E);
        variable.Eang = angle(E);
        variable.Iabs = abs(I);
        variable.I = I;


% Get data from Interfacing Converters
% [variable.IC_P,variable.IC_Q,variable.IC_Edc] = Get_IC_data(Grid_para);

end