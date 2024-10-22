function [] =  sendToZenone(setpoint,Grid_para)

% setpoint.P = 1000
% setpoint.Q = 1000

port = 35011;
host = '192.168.1.12';

sp_P = setpoint.P*Grid_para.A_b; 
sp_Q = setpoint.Q*Grid_para.A_b; 

data = struct('P', sp_P, 'Q', sp_Q);

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
% 
% 