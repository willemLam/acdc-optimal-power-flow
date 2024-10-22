function [SM_info,variables] = ReceiveFromCC(Grid_para,variables,Topology) 

if Topology == "ISLAND"
    port = 35002;
    host = '192.168.1.210';
elseif Topology == "CONNECTED"
    port = 35081;
    host = '192.168.1.210';
else
    error()
end

count = 0;

try
    [message, ~] = judp('RECEIVE',port,2^15,4000); 
    json1 = char(message);
    all = loadjson(json1');

    SM_info.mode = all.Data.Operation_mode;
    SM_info.Optimal_slack = all.Data.Optimal_slack;
    SM_info.Angle = all.Data.upper_layer_grid_info.Angle/(180/pi);
    SM_info.Frequency = all.Data.upper_layer_grid_info.Frequency;
    SM_info.Magnitude = all.Data.upper_layer_grid_info.Magnitude/(Grid_para.V_b/sqrt(3));
    SM_info.SC_soc = all.Data.SC_soc;

%% from MeteoBox
    variables.PVmax_P_perun = 0*max(0,all.Data.PV.perun)/Grid_para.A_b;  
    variables.PVmax_P_solarmax = max(0,all.Data.PV.solarmax)/Grid_para.A_b; 
    variables.PVmax_P_facade = max(0,all.Data.PV.facade)/Grid_para.A_b;  
%     variables.PVmax_P_emul4 = 0.03;%max(0,all.Data.PV.emul4)/Grid_para.A_b;
%     variables.PVmax_P_emul5 = 0.03;%max(0,all.Data.PV.emul5)/Grid_para.A_b;
%     variables.PVmax_P_emul6 = 0.03;%max(0,all.Data.PV.emul6)/Grid_para.A_b;
%     variables.PVmax_P_emul7 = 0.03;%max(0,all.Data.PV.emul7)/Grid_para.A_b;
catch
    count = count + 1;
    if count == 3
           warning('Cannot connect to State Machine')
    end
end




end