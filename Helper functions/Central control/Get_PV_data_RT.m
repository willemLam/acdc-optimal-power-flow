function PVmax = Get_PV_data_RT(TYPE,Grid_para,P_perun_prev, P_solarmax_prev, P_facade_prev, P_PV2_prev, P_PV3_prev)


count = 0;
try
    if TYPE == 'EXPERIMENT' %'SIMULATION'
        UTC = 2;
        tf = 1694269322 + 10000; %  matlabToEpochTime(now())-UTC*3600 -1e8 %- 3600*100 + (i-1)*5*60; 
        t0 = tf-2*3600*24;
        [P_max_perun, P_max_solarmax, P_max_facade] = PV_history(tf,t0);


    elseif TYPE == 'SIMULATION'
        P_max_perun = 0;
        P_max_solarmax = 0;
        P_max_facade = 0 ;
        %[P_max_perun, P_max_solarmax, P_max_facade] = readFromMeteoBox();
    end
    PVmax.P_perun = P_max_perun/Grid_para.A_b;  % Perun is off
    PVmax.P_solarmax = P_max_solarmax/Grid_para.A_b;
    PVmax.P_facade = P_max_facade/Grid_para.A_b;
    PVmax.P_emul4 = 0;
    PVmax.P_emul5 = 0;
    PVmax.P_emul6 = 0;
    PVmax.P_emul7 = 0;

catch
    disp('Cannot connect to Influx')
    count = count + 1;
    if count == 2
        PVmax.P_perun = 0;  %Perun is off
        PVmax.P_solarmax = 0; %previous PV (accounting for the curtailment)
        PVmax.P_facade = 0; %previous PV 
        PVmax.P_emul4 = 0;
        PVmax.P_emul5 = 0;
        PVmax.P_emul6 = 0;
        PVmax.P_emul7 = 0;
        disp('Previous GHI is used')
    end
end

end
