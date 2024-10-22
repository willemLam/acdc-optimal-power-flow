function dF = Mismatch_1ph_V2_quadratic_loss(E,S,E_star,S_star,Grid_para,Filter_para,idx)
    
    n_dc = Grid_para.n_dc;
    n_ac = Grid_para.n_ac;
    n_nodes = Grid_para.n_nodes;
    n_ph = Grid_para.n_ph;
    G = Grid_para.G;
    B = Grid_para.B;
    Y = complex(G,B);    
    I_b =Grid_para.Y_b*Grid_para.V_b;
    
    pos_ac3 = Grid_para.pos_ac3;
    pos_dc3 = Grid_para.pos_dc3;
    

    %% compute quadratic losses
    % Pac + Ploss = Pdc
    Iac = abs(Y(pos_ac3(:,1),:) * E);
    Idc = abs(Y(pos_dc3(:,1),:) * E);
    Vdc = E(pos_dc3(:,1));
    P_quad_ic =                  Filter_para.b .* Iac + Filter_para.c .* Iac.^2 ...
               + Filter_para.a .*Vdc + Filter_para.d .* Idc + Filter_para.e .* Idc.^2;


    if Filter_para.Include_losses == 0
        P_quad_ic = zeros(length(pos_ac3(:,1)),1);
    end
    

    P_quad_ic_PQ = P_quad_ic(Grid_para.AFE_type=="PQ");
    P_quad_ic_VQ = P_quad_ic(Grid_para.AFE_type=="VQ");
    P_quad_ic_VV = P_quad_ic(Grid_para.AFE_type=="VV");

    %% Vdc - Q nodes
    % Recompute the DC voltage at the AFE


    for t = 1:length(idx.vscdc_vq)
        
        p = find(pos_dc3(:,1)==idx.vscdc_vq(t));
        alphan =  ( G(pos_dc3(p,1),pos_dc3(p,2))* E(pos_dc3(p,2)) )^2 - 4 * G(pos_dc3(p,1),pos_dc3(p,1)) * real( S(idx.vscac_vq(p)) + P_quad_ic(p) );
        E(idx.vscdc_vq(p)) = -G(pos_dc3(p,1),pos_dc3(p,2))/(2*G(pos_dc3(p,1),pos_dc3(p,1))).* E(pos_dc3(p,2)) + sqrt(alphan)/(2*G(pos_dc3(p,1),pos_dc3(p,1)));
    end
        

    %% Vac - theta nodes
    for t = 1:length(idx.vscdc_vv)
        
        p = find(pos_dc3(:,1)==idx.vscdc_vv(t));
        
        %solve to E_l real
        a = G(pos_ac3(p,1),pos_ac3(p,1));
        b = real(conj(Y(pos_ac3(p,1),pos_ac3(p,2))* (E(pos_ac3(p,2)))));
        c = + G(pos_ac3(p,1),pos_ac3(p,1)) * imag(E(pos_ac3(p,1))).^2  ...
            - imag(E(pos_ac3(p,1))) .* imag(conj(Y(pos_ac3(p,1),pos_ac3(p,2))* E(pos_ac3(p,2)))) ...
            + real(S(pos_dc3(p,1)) + P_quad_ic(p)) ; %P_quad_ic is dependent on Vac, but we ignore it
    
        % check_ should be 0
        % a * real(E(pos_ac3(p,1))).^2 + b .* real(E(pos_ac3(p,1))) + c
        
        %E_l real
        D = sqrt(b.^2-4 .*a .*c);
        
        E_l_real = real((-b + D ) ./ (2*a));
        
        %solve to E_l imag
        a = diag(G(pos_ac3(p,1),pos_ac3(p,1)));
        b = -imag(conj(Y(pos_ac3(p,1),pos_ac3(p,2))* (E(pos_ac3(p,2)))));
        c = + G(pos_ac3(p,1),pos_ac3(p,1)) * real(E(pos_ac3(p,1))).^2 ...
            + real(E(pos_ac3(p,1))) .* real(conj(Y(pos_ac3(p,1),pos_ac3(p,2))* (E(pos_ac3(p,2))))) ...
            + real(S(pos_dc3(p,1)) - P_quad_ic(p));
        % check_ should be 0   
        % a .* imag(E(pos_ac3(p,1))).^2 + b .* imag(E(pos_ac3(p,1))) + c
        
        %E_l imag
%         imag(E(pos_ac3(p,1)))
        D = sqrt(b.^2-4 .*a .*c);
        E_l_imag = real((-b + D ) ./ (2*a));
        
        E(idx.vscac_vv(t)) = complex(E_l_real,E_l_imag);
       
    end
    
    %% P - Q nodes
    P_ic = real(S_star);
    for t = 1:length(idx.vscac_pq)   
        p = find(pos_ac3(:,1)==idx.vscac_pq(t));
        P_ic(pos_ac3(p,1)) = P_ic(pos_ac3(p,1)); % Pdc is the reference
        P_ic(pos_dc3(p,1)) = -P_ic(pos_ac3(p,1)) - P_quad_ic(p); %Pdc Pac = -Pdc - Ploss
        
    end
    S_star_ic = complex(P_ic,imag(S_star));

    % Compute the mismatches for the entire network.
    dS = S_star_ic-S; %include losses
    dP = real(dS);
    dQ = imag(dS);
    dV2 = abs(E_star).^2-abs(E).^2; 
    dE = E_star - E; 
    dPdc = real(dS);
    dEdc = real(E_star-E);

    
    
    dF = [ dP(idx.pqac); 
           dQ(idx.pqac); 
           dP(idx.pvac); 
           dV2(idx.pvac);              
           dP(idx.vscac_pq);
           dQ(idx.vscac_pq);
           dEdc(idx.vscdc_vq);
           dQ(idx.vscac_vq);
           real(dE(idx.vscac_vv));
%            imag(dE(idx.vscac_vv)); %not needed
           dPdc(idx.vscdc_pq);
           dPdc(idx.pdc)];
     
end