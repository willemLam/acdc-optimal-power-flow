function [] =  sendToZenone(P, Q, time, iter, dest_port, dest_IP)
% write to a specif host and port, with the
% Qac, Vdc in (kVar, V)
% defined ports and IP for the resources 



%--- example of arguments
% P = 1000
% Q = 1000
% time = 6
% iter = 1
% port = 35011;
% host = '192.168.1.12';

port = dest_port;
host = dest_IP;

u = udp(host,port);
fopen(u);

for k = 1:3 % send mulitple times to make sure UDP packet is recieved.
    
    set_points = struct('P', P, 'Q', Q);
    % set_points = struct('Vdc1', Vdc1, 'Qac1', Qac1,'iteration', iter, 'time', time);
    json_str = savejson('',set_points');

    fwrite(u, json_str)
%     pause(1)

end

fclose(u);

end
% 
% 