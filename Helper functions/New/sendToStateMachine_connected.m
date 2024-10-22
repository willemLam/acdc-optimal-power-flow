function SM_info = sendToStateMachine_connected(SM_info)

port = 35234;
host = '192.168.1.210';

count = 0;

    SM_info.Operation_mode = 'Grid connected';
    SM_info.Optimal_slack = 'AFE 1';

Info = struct('Resource', 'Connected OPF', 'Problem', 0, 'Timestamp',SM_info.timestep);
grid = struct('Frequency',0, 'Angle', 0, 'Magnitude', 0);
Data = struct('Operation_mode', SM_info.Operation_mode, 'solution_feasible', SM_info.feasibility, 'Optimal_slack', SM_info.Optimal_slack, 'setpoints', SM_info.setpoints, 'upper_layer_grid_info', grid);

all = struct('Info', Info,  'Data', Data);

json_str = savejson('',all');

u = udp(host,port,OutputBufferSize=2^16);
fopen(u);


try
    fwrite(u, json_str)
catch
    count = count + 1;
    if count == 3
           warning('Cannot connect to State Machine')
    end
end


fclose(u);


end