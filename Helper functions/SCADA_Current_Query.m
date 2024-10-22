function Current = SCADA_Current_Query(Current_idx,t0, tf, resolution,Grid_para)

    ac_nodes = 1:Grid_para.n_ac;
    dc_nodes = Grid_para.n_ac+1:Grid_para.n_ac + Grid_para.n_dc;

    
    for i = 1:length(Current_idx)

        if ismember(Current_idx(i),ac_nodes)
            field = string(sprintf('B%02dIA_mag',Current_idx(i)));
            [tab_A] = fluxQuery('microgrid_ST', t0, tf, 'microgrid', 'Resource', 'SCADA', field, resolution);

            field = string(sprintf('B%02dIB_mag',Current_idx(i)));
            [tab_B] = fluxQuery('microgrid_ST', t0, tf, 'microgrid', 'Resource', 'SCADA', field, resolution);

            field = string(sprintf('B%02dIC_mag',Current_idx(i)));
            [tab_C] = fluxQuery('microgrid_ST', t0, tf, 'microgrid', 'Resource', 'SCADA', field, resolution);

            newTimetable_Imag = synchronize(tab_A, tab_B, tab_C);

            
 
        elseif ismember(Current_idx(i),dc_nodes)
            field = string(sprintf('B%02dIA_mag',Current_idx(i)));
            [newTimetable_Imag] = fluxQuery('microgrid_ST', t0, tf, 'microgrid', 'Resource', 'SCADA', field, resolution);
            
        else
            warning('wrong')
        end

        if i == 1
            Current = newTimetable_Imag;
        else
            Current = synchronize(Current,newTimetable_Imag);
        end

    end

end