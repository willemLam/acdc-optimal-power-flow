function status_switch = Read_from_islander()

    UDP_port = 35692;
    [message, ~] = judp('RECEIVE',UDP_port,2^18); 
    message = loadjson(char(message));
    status_switch = message.data.status_island;
end
