function [P_max_perun, P_max_solarmax, P_max_facade, P_max_PV2, P_max_PV3] = Get_PV_data(TYPE,Grid_para,P_perun_prev, P_solarmax_prev, P_facade_prev, P_PV2_prev, P_PV3_prev)

A_b = Grid_para.A_b;

count = 0;
try
    if TYPE == 'SIMULATION'
        UTC = 2;
        tf = 1694269322; %matlabToEpochTime(now)-UTC*3600 - 3600*100 + (i-1)*5*60;
        t0 = tf-2*3600*24;
        [P_max_perun, P_max_solarmax, P_max_facade] = PV_history(tf,t0);

    elseif TYPE == 'EXPERIMENT'
        [P_max_perun, P_max_solarmax, P_max_facade] = PV_forecast_update();
    end
    P_max_perun = P_max_perun/A_b;  % Perun is off
    P_max_solarmax = P_max_solarmax/A_b;
    P_max_facade = P_max_facade/A_b;
    P_max_PV2 = 0;
    P_max_PV3 = 0;

catch
    disp('Cannot connect to Influx')
    count = count + 1;
    if count == 2
        P_max_perun = double(P_perun_prev);  %Perun is off
        P_max_solarmax = double(P_solarmax_prev); %previous PV (accounting for the curtailment)
        P_max_facade = double(P_facade_prev); %previous PV 
        disp('Previous GHI is used')
    end
end

end
