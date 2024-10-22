function [Kp_abs, Kq_abs, Kvabs_abs,Kvang_abs, Kp_ang, Kq_ang, Kvabs_ang, Kvang_ang] = transform_K_polar(K,Grid_para, idx)

    J_P_abs = zeros(Grid_para.n_nodes);
    J_Q_abs = zeros(Grid_para.n_nodes);
    J_Vabs_abs = zeros(Grid_para.n_nodes);
    J_Vang_abs = zeros(Grid_para.n_nodes);
    J_P_ang = zeros(Grid_para.n_nodes);
    J_Q_ang = zeros(Grid_para.n_nodes);
    J_Vabs_ang = zeros(Grid_para.n_nodes);
    J_Vang_ang = zeros(Grid_para.n_nodes);
    for k = 1:size(K,1)
        if( sum( K{k,1} == idx.slack))
            continue
        elseif( sum( K{k,1} == idx.pqac))
            J_P_abs(:,k) = (K{k,5}{1,1});
            J_Q_abs(:,k) = (K{k,5}{2,1});
            J_P_ang(:,k) = (K{k,6}{1,1});
            J_Q_ang(:,k) = (K{k,6}{2,1});

        elseif( sum( K{k,1} == idx.pvac ))
            J_P_abs(:,k) =    (K{k,5}{1,1});
            J_Vabs_abs(:,k) = (K{k,5}{2,1});
            J_P_ang(:,k) =    (K{k,6}{1,1});
            J_Vabs_ang(:,k) = (K{k,6}{2,1});

        elseif( sum( K{k,1} == idx.pdc ) )
            J_P_abs(:,k) = (K{k,5}{1,1});
            J_P_ang(:,k) = (K{k,6}{1,1});

        elseif( sum( K{k,1} == idx.vdc ) )
            J_Vabs_abs(:,k) = (K{k,5}{1,1});
            J_Vabs_ang(:,k) = (K{k,6}{1,1});

        elseif( sum( K{k,1} == idx.vscac_pq))
            J_P_abs(:,k) = (K{k,5}{1,1});    
            J_Q_abs(:,k) = (K{k,5}{2,1});
            J_P_ang(:,k) = (K{k,6}{1,1});    
            J_Q_ang(:,k) = (K{k,6}{2,1});

        elseif( sum( K{k,1} == idx.vscac_vq))
            J_P_abs(:,k) = (K{k,5}{1,1});
            J_Q_abs(:,k) = (K{k,5}{2,1});
            J_P_ang(:,k) = (K{k,6}{1,1});
            J_Q_ang(:,k) = (K{k,6}{2,1});

        elseif( sum( K{k,1} == idx.vscac_vv))
            J_Vabs_abs(:,k) = (K{k,5}{1,1});
            J_Vang_abs(:,k) = (K{k,5}{2,1});
            J_Vabs_ang(:,k) = (K{k,6}{1,1});
            J_Vang_ang(:,k) = (K{k,6}{2,1});

        elseif( sum( K{k,1} == idx.vscdc_pq ))
            J_P_abs(:,k) = (K{k,5}{1,1});
            J_P_ang(:,k) = (K{k,6}{1,1});

        elseif( sum( K{k,1} == idx.vscdc_vq ))
            J_Vabs_abs(:,k) = (K{k,5}{1,1});
            J_Vabs_ang(:,k) = (K{k,6}{1,1});

        elseif( sum( K{k,1} == idx.vscdc_vv ))
            J_Vabs_abs(:,k) = (K{k,5}{1,1});
            J_Vabs_ang(:,k) = (K{k,6}{1,1});

        else
            warning('somethings off mate')
        end
    end

    Kp_abs = J_P_abs;%(4:end,4:end);
    Kq_abs = J_Q_abs;%(4:end,4:end);
    Kvabs_abs = J_Vabs_abs;%(4:end,4:end);
    Kvang_abs = J_Vang_abs;%(4:end,4:end);
    Kp_ang = J_P_ang;%(4:end,4:end);
    Kq_ang = J_Q_ang;%(4:end,4:end);
    Kvabs_ang = J_Vabs_ang;%(4:end,4:end);
    Kvang_ang = J_Vang_ang;%(4:end,4:end);
end