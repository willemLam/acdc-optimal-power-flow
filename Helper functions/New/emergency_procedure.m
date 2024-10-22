function emergency_procedure(Grid_para)

%% AFE
        count = 0;

        setpoint.Q = [0; 0; 0; 0];
        setpoint.Edc = [1; 1; 1; 1];
        setpoint.Eac = [1; 1];
        setpoint.freq = [50; 50];
        try
            sendToAFE(setpoint,Grid_para)
        catch
            count = count + 1;
            if count == 3
                   warning('Cannot connect to AFE')
            end
        end

%% PV emulator and SuperCap
        count = 0;

        setpoint.P = [0; 0; 0; 0];
        try
            sendToPV_SC(setpoint,Grid_para)
        catch
            count = count + 1;
            if count == 3
                   warning('Cannot connect to PV emulator and Super Cap')
            end
        end
%% Zenone
        count = 0;

        setpoint.P = 0;
        setpoint.Q = 0;
        try
            sendToZenone(setpoint,Grid_para)
        catch
            count = count + 1;
            if count == 3
                   warning('Cannot connect to Zenone')
            end
        end

end
