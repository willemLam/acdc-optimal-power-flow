function Syncro = SCADA_Syncro_Query(Syncro_idx,t0, tf, resolution,Grid_para)

    alp = exp(2*pi/3*1i);
    T_inv = 1/3*   [1 1     1; 
                    1 alp   alp^2; 
                    1 alp^2 alp];
    fields = string({'DOWN_VA_mag'; 'DOWN_VB_mag'; 'DOWN_VC_mag';...
             'UP_VA_mag';'UP_VB_mag';'UP_VC_mag';...
             'DOWN_VA_ph';'DOWN_VB_ph';'DOWN_VC_ph';...
             'UP_VA_ph';'UP_VB_ph';'UP_VC_ph';...
             'DOWN_IA_mag';'DOWN_IB_mag';'DOWN_IC_mag';...
             'f_downstream';'f_upstream'});

    for i = 1:length(fields)

       [newTimetable] = fluxQuery('microgrid_ST', t0, tf, 'microgrid', 'Resource', 'SCADA', fields(i), resolution);

        if i == 1
            Syncro = newTimetable;
        else
            Syncro = synchronize(Syncro,newTimetable);
        end

    end

    %make symetric components

    direct_component  = [Syncro.UP_VA_mag Syncro.UP_VB_mag Syncro.UP_VC_mag] .* exp(1i*[Syncro.UP_VA_ph Syncro.UP_VB_ph Syncro.UP_VC_ph]) * T_inv(:,2);
    Syncro.UP_V_mag = abs(direct_component);
    Syncro.UP_V_ph  = angle(direct_component);

    direct_component  = [Syncro.DOWN_VA_mag Syncro.DOWN_VB_mag Syncro.DOWN_VC_mag] .* exp(1i*[Syncro.DOWN_VA_ph Syncro.DOWN_VB_ph Syncro.DOWN_VC_ph]) * T_inv(:,2);
    Syncro.DOWN_V_mag = abs(direct_component);
    Syncro.DOWN_V_ph  = angle(direct_component);

end