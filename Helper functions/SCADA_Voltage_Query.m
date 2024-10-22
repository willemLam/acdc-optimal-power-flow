function Voltage = SCADA_Voltage_Query(Voltage_idx,t0, tf, resolution,Grid_para)

    ac_nodes = 1:22;%1:Grid_para.n_ac;
    dc_nodes = 23:30;%Grid_para.n_ac+1:Grid_para.n_ac + Grid_para.n_dc;

    alp = exp(2*pi/3*1i);
    T_inv = 1/3*   [1 1     1; 
                    1 alp   alp^2; 
                    1 alp^2 alp];

    for i = 1:length(Voltage_idx)

        if ismember(Voltage_idx(i),ac_nodes)
            field = string(sprintf('B%02dVA_mag',Voltage_idx(i)));
            [tab_A] = fluxQuery('microgrid_ST', t0, tf, 'microgrid', 'Resource', 'SCADA', field, resolution);

            field = string(sprintf('B%02dVB_mag',Voltage_idx(i)));
            [tab_B] = fluxQuery('microgrid_ST', t0, tf, 'microgrid', 'Resource', 'SCADA', field, resolution);

            field = string(sprintf('B%02dVC_mag',Voltage_idx(i)));
            [tab_C] = fluxQuery('microgrid_ST', t0, tf, 'microgrid', 'Resource', 'SCADA', field, resolution);

            alignedTimetables_mag = synchronize(tab_A, tab_B, tab_C);
            
            field = string(sprintf('B%02dVA_ph',Voltage_idx(i)));
            [tab_A] = fluxQuery('microgrid_ST', t0, tf, 'microgrid', 'Resource', 'SCADA', field, resolution);

            field = string(sprintf('B%02dVB_ph',Voltage_idx(i)));
            [tab_B] = fluxQuery('microgrid_ST', t0, tf, 'microgrid', 'Resource', 'SCADA', field, resolution);

            field = string(sprintf('B%02dVC_ph',Voltage_idx(i)));
            [tab_C] = fluxQuery('microgrid_ST', t0, tf, 'microgrid', 'Resource', 'SCADA', field, resolution);

            alignedTimetables_ph = synchronize(tab_A, tab_B, tab_C);

            direct_component  = alignedTimetables_mag.Variables .* exp(1i*alignedTimetables_ph.Variables) * T_inv(:,2);
            
            
            newTimetable = timetable(alignedTimetables_mag.time, abs(direct_component), angle(direct_component) , 'VariableNames', {sprintf('B%02dV_abs',Voltage_idx(i)),sprintf('B%02dV_ph',Voltage_idx(i))});
    

        elseif ismember(Voltage_idx(i),dc_nodes)
            field = string(sprintf('B%02dVA_mag',Voltage_idx(i)));
            [tab] = fluxQuery('microgrid_ST', t0, tf, 'microgrid', 'Resource', 'SCADA', field, resolution);
            newTimetable = timetable(tab.time, tab.Variables, 'VariableNames', string(sprintf('B%02dV_mag',Voltage_idx(i))));
    
        else
            warning('wrong')
        end

        if i == 1
            Voltage = newTimetable;
        else
            Voltage = synchronize(Voltage,newTimetable);
        end

    end

end