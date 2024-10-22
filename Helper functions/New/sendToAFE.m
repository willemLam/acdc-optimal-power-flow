function sendToAFE(setpoint,Grid_para)

% setpoint.Q = [0.05; 0.05; 0.05; 0.05];
% setpoint.Edc = [0.9; 0.95; 1; 1.05];
% setpoint.Eac = [0.9+0.1*1i; 0.95+0.5*1i];
% setpoint.freq = [49; 51];

port = 34720;
host = "192.168.1.20";

count = 0;

sp_Q = setpoint.Q*Grid_para.A_b; 
sp_Q(abs(sp_Q) < 100) = 0; % no need to inject values <100
sp_Edc = setpoint.Edc*Grid_para.Vdc_b;
sp_Eac = setpoint.Eac*Grid_para.V_b/sqrt(3);
sp_freq = setpoint.freq;

data.references = struct('Vdc1', sp_Edc(1), 'Qac1', sp_Q(1), 'V_abs1', abs(sp_Eac(1)), 'V_arg1', angle(sp_Eac(1)), 'freq1', sp_freq(1), 'Vdc2', sp_Edc(2), 'Qac2', sp_Q(2), 'V_abs2', abs(sp_Eac(2)), 'V_arg2', angle(sp_Eac(2)), 'freq2', sp_freq(2),'Vdc3', sp_Edc(3), 'Qac3', sp_Q(3),'Vdc4', sp_Edc(4), 'Qac4', sp_Q(4));
data.Start_forming_1 = true;
data.Start_forming_2 = false;

json_str = savejson('',data');

u = udp(host,port,OutputBufferSize=2^16);
fopen(u);


try
    fwrite(u, json_str)
catch
    count = count + 1;
    if count == 3
           warning('Cannot connect to Interfacing converters')
    end
end


fclose(u);

end
