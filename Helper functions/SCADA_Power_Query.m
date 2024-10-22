function Power = SCADA_Power_Query(Power_idx,t0, tf, resolution,Grid_para)

    ac_nodes = 1:Grid_para.n_ac;
    dc_nodes = Grid_para.n_ac+1:Grid_para.n_ac + Grid_para.n_dc;

    for i = 1:length(Power_idx)

        if ismember(Power_idx(i),ac_nodes)
            field = string(sprintf('B%02dPA',Power_idx(i)));
            [tab_A] = fluxQuery('microgrid_ST', t0, tf, 'microgrid', 'Resource', 'SCADA', field, resolution);

            field = string(sprintf('B%02dPB',Power_idx(i)));
            [tab_B] = fluxQuery('microgrid_ST', t0, tf, 'microgrid', 'Resource', 'SCADA', field, resolution);

            field = string(sprintf('B%02dPC',Power_idx(i)));
            [tab_C] = fluxQuery('microgrid_ST', t0, tf, 'microgrid', 'Resource', 'SCADA', field, resolution);

            alignedTimetables = synchronize(tab_A, tab_B, tab_C);
            summedValues = sum(alignedTimetables.Variables, 2);
    
            newTimetable_P = timetable(alignedTimetables.time, summedValues, 'VariableNames', string(sprintf('B%02dP',Power_idx(i))));
    
            field = string(sprintf('B%02dQA',Power_idx(i)));
            [tab_A] = fluxQuery('microgrid_ST', t0, tf, 'microgrid', 'Resource', 'SCADA', field, resolution);

            field = string(sprintf('B%02dQB',Power_idx(i)));
            [tab_B] = fluxQuery('microgrid_ST', t0, tf, 'microgrid', 'Resource', 'SCADA', field, resolution);

            field = string(sprintf('B%02dQC',Power_idx(i)));
            [tab_C] = fluxQuery('microgrid_ST', t0, tf, 'microgrid', 'Resource', 'SCADA', field, resolution);

            alignedTimetables = synchronize(tab_A, tab_B, tab_C);
            summedValues = sum(alignedTimetables.Variables, 2);
    
            newTimetable_Q = timetable(alignedTimetables.time, summedValues, 'VariableNames', string(sprintf('B%02dQ',Power_idx(i))));
    
            newTimetable = synchronize(newTimetable_P,newTimetable_Q);

        elseif ismember(Power_idx(i),dc_nodes)
            field = string(sprintf('B%02dPA',Power_idx(i)));
            [tab] = fluxQuery('microgrid_ST', t0, tf, 'microgrid', 'Resource', 'SCADA', field, resolution);
            newTimetable = timetable(tab.time, tab.Variables, 'VariableNames', string(sprintf('B%02dP',Power_idx(i))));
    
        else
            warning('wrong')
        end

        if i == 1
            Power = newTimetable;
        else
            Power = synchronize(Power,newTimetable);
        end

    end

end