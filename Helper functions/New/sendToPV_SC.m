function sendToPV_SC(setpoint,Grid_para)

% setpoint.P = [-0.005; 0.005; -0.005; 0.005];


port = 37690;
host = "192.168.1.73";

count = 0;

sp_P = setpoint.P*Grid_para.A_b; 

data = struct('P1', sp_P(1), 'P2', sp_P(2), 'P3', sp_P(3), 'P4', sp_P(4));

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
