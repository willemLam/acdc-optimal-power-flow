function variables = Receive_data_from_CC(Grid_para,variables)

% info


CC_port = 35211;
CC_host = '192.168.1.86';

SCADA_port = 35005;
SCADA_host = '192.168.1.2';

count = 0;

try

[message, ~] = judp('RECEIVE',CC_port,2^15,4000); 
json1 = char(message);
data = loadjson(json1');
% data = jsondecode(json1');

variables.PV_P_facade = data.Facade.Data.P /Grid_para.A_b;

catch
    count = count + 1;
    if count == 3
           warning('Cannot connect to CC')
    end
end


end