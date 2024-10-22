function [data] = readFromResources(reading_port, byte_size )

[message, ~] = judp('RECEIVE', reading_port, byte_size);
json1 = char(message);
% fprintf(1,'%s\n',json1);
data=loadjson(json1');

end
