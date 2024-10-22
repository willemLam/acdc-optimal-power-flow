function variables = Receive_raw_data_from_SCADA(Grid_para,variables,prev_solution)

% info


CDA_port = 35210;
CDA_host = '192.168.1.86';

SCADA_port = 35005;
SCADA_host = '192.168.1.2';

count = 0;

try


    [message, ~] = judp('RECEIVE',35005,2^18,4000); 
    json1 = char(message);
%     fprintf(1,'%s\n',json1);
    data = loadjson(json1');

    measurements = data.Data.measurements;

    variables.Eac_raw = complex(measurements(1:21),measurements(22:42))*Grid_para.V_b;
    variables.Iac_raw = complex(measurements(43:54),measurements(55:66))*Grid_para.I_b;
    variables.Edc_raw = measurements(247:254)*Grid_para.Vdc_b;
    variables.Idc_raw = measurements(255:262)*Grid_para.Idc_b;


catch
    count = count + 1;
    if count == 3
           warning('Cannot connect to SCADA')
    end

    variables = prev_solution;
end


end