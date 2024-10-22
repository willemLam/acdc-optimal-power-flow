function Tsolution_combined = Get_all_data_old(FileName)
    
    s = whos('-file',FileName);
    N = length({s.name});
    solution = struct();
    
    for i = 1:N
      
            if rem(i,100) == 0 %update on loading the data
                i
            end
        name = ['data_',num2str(i)];
        D = load(FileName, genvarname(name));
        d = D.(genvarname(name));
    
        solution.PVmax_P_perun(:,i) = d.PVmax_P_perun;
        solution.PVmax_P_solarmax(:,i) = d.PVmax_P_solarmax;
        solution.PVmax_P_facade(:,i) = d.PVmax_P_facade;
    
        solution.PV_P_perun(:,i) = d.PV_P_perun;
        solution.PV_P_solarmax(:,i) = d.PV_P_solarmax;
        solution.PV_P_facade(:,i) = d.PV_P_facade;
    
        solution.P(:,i) = d.P;
        solution.Q(:,i) = d.Q;
        solution.S(:,i) = d.S;
        solution.E(:,i) = d.E;
        solution.Iabs(:,i) = d.Iabs;
    
        solution.reconstructed.S(:,i) = d.reconstructed.S;
        solution.reconstructed.E(:,i) = d.reconstructed.E;
        solution.reconstructed.Iabs(:,i) = d.reconstructed.Iabs;
    
        solution.Ploss(:,i) = d.Ploss;
        solution.Qloss(:,i) = d.Qloss;
        solution.Sambat_P(:,i) = d.Sambat_P;
        solution.SupCap_P(:,i) = d.SupCap_P;
        solution.SupCap_SOC(:,i) = d.SupCap_SOC*100;
        solution.IC_Q(:,i) = d.IC_Q;
    %     solution.IC_E(:,i) = d.IC_E;
        solution.feasible(i) = d.proceed;
    
        solution.mode{:,i} = d.mode;
        solution.stime{:,i} = d.time;
    
        solution.duration(:,i) = d.duration;
    
        variables.E(:,i) = d.variables.E;
        variables.S(:,i) = d.variables.S;
        variables.I(:,i) = d.variables.I;
        solution.f_downstream(:,i) = d.variables.f_downstream;
        solution.f_upstream(:,i) = d.variables.f_upstream;
        solution.V_downstream(:,i) = d.variables.V_downstream;
        solution.V_upstream(:,i) = d.variables.V_upstream;

    end

        % Start with the time array
%     time_array = datetime([solution.stime{:}])';
%     
%     % Initialize the combined timetable with time array
%     Tsolution_combined = timetable(time_array, 'RowTimes', time_array);

    % List of fields to add to the combined timetable
    fields_to_convert = {'PVmax_P_perun', 'PVmax_P_solarmax', 'PVmax_P_facade', ...
                         'PV_P_perun', 'PV_P_solarmax', 'PV_P_facade', ...
                         'P', 'Q', 'S', 'E', 'Iabs', ...
                         'Ploss', 'Qloss', 'Sambat_P', 'SupCap_P', 'SupCap_SOC', 'IC_Q', ...
                         'duration','f_downstream','f_upstream','V_downstream','V_upstream'};
    
    % Add each field to the combined timetable
    for k = 1:numel(fields_to_convert)
        field_name = fields_to_convert{k};
        Tsolution_combined.(field_name) = solution.(field_name)';
    end
    
    % Add reconstructed fields
    reconstructed_fields = {'S', 'E', 'Iabs'};
    for k = 1:numel(reconstructed_fields)
        field_name = ['reconstructed_' reconstructed_fields{k}];
        Tsolution_combined.(field_name) = solution.reconstructed.(reconstructed_fields{k})';
    end
    
    % Add the feasible and mode fields
    Tsolution_combined.feasible = solution.feasible';
    Tsolution_combined.mode = solution.mode';
    
    % Add variables fields
    variables_fields = {'E', 'S', 'I'};
    for k = 1:numel(variables_fields)
        field_name = ['variables_' variables_fields{k}];
        Tsolution_combined.(field_name) = variables.(variables_fields{k})';
    end


end