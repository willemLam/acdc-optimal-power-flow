function [Kp, Kq, Kvabs, Kvang] = transform_K_complex(K,Grid_para, idx)

    J_P = zeros(Grid_para.n_nodes);
    J_Q = zeros(Grid_para.n_nodes);
    J_Vabs = zeros(Grid_para.n_nodes);
    J_Vang = zeros(Grid_para.n_nodes);

    for k = 1:size(K,1)
        if( sum( K{k,1} == idx.slack))
            continue
        elseif( sum( K{k,1} == idx.pqac))
            J_P(:,k) = (K{k,2}{1,1});
            J_Q(:,k) = (K{k,2}{2,1});

        elseif( sum( K{k,1} == idx.pvac ))
            J_P(:,k)    = (K{k,2}{1,1});
            J_Vabs(:,k) = (K{k,2}{2,1});

        elseif( sum( K{k,1} == idx.pdc ) )
            J_P(:,k) = (K{k,2}{1,1});
         
        elseif( sum( K{k,1} == idx.vdc ) )
            J_Vabs(:,k) = (K{k,2}{1,1});

        elseif( sum( K{k,1} == idx.vscac_pq))
            J_P(:,k) = (K{k,2}{1,1});    
            J_Q(:,k) = (K{k,2}{2,1});

        elseif( sum( K{k,1} == idx.vscac_vq))
            J_P(:,k) = (K{k,2}{1,1});
            J_Q(:,k) = (K{k,2}{2,1});

        elseif( sum( K{k,1} == idx.vscac_vv))
            J_Vabs(:,k) = (K{k,2}{1,1});
            J_Vang(:,k) = (K{k,2}{2,1});

        elseif( sum( K{k,1} == idx.vscdc_pq ))
            J_P(:,k) = (K{k,2}{1,1});

        elseif( sum( K{k,1} == idx.vscdc_vq ))
            J_Vabs(:,k) = (K{k,2}{1,1});
        
        elseif( sum( K{k,1} == idx.vscdc_vv ))
            J_P(:,k) = (K{k,2}{1,1});
        else
            warning('somethings off mate')
        end
    end

    Kp = J_P;%(4:end,4:end);
    Kq = J_Q;%(4:end,4:end);
    Kvabs = J_Vabs;%(4:end,4:end);
    Kvang = J_Vang;
end