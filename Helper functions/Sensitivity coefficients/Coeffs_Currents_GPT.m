function [K_Iabs_P, K_I_P, K_Iabs_Q, K_I_Q, K_Iabs_Eabs, K_I_Eabs, K_Iabs_Eang, K_I_Eang ] = Coeffs_Currents_GPT(I, E, K_E, Grid_para)
    YYL = Grid_para.YYL;
    YYT = Grid_para.YYT;
    nph = Grid_para.n_ph;
   
    list_i = Grid_para.lines(:,1);
    list_j = Grid_para.lines(:,2);

    num_power_injections = size(K_E.P, 2);
    num_lines = length(list_i);

    K_Iabs_P = zeros(num_lines, num_power_injections);
    K_Iabs_Q = zeros(num_lines, num_power_injections);
    K_Iabs_Eabs = zeros(num_lines, num_power_injections);
    K_Iabs_Eang = zeros(num_lines, num_power_injections);

    K_I_P = zeros(num_lines, num_power_injections);
    K_I_Q = zeros(num_lines, num_power_injections);
    K_I_Eabs = zeros(num_lines, num_power_injections);
    K_I_Eang = zeros(num_lines, num_power_injections);
   
    tol = 1e-8;

    % Compute IKp_abs, IKq_abs, and IKv_abs
    for n = 1:num_power_injections
        n3 = (n-1)*nph+1 : n*nph;

        for l = 1:num_lines
            l3 = (l-1)*nph+1 : l*nph;
            i = (list_i(l)-1)*nph+1 : list_i(l)*nph;
            j = (list_j(l)-1)*nph+1 : list_j(l)*nph;

            Delta_Eabs = abs(abs(E(list_i(l))) - abs(E(list_j(l))));
            Delta_Eang = abs(angle(E(list_i(l))) - angle(E(list_j(l))));
            
            if any(abs(I(l3))) &&  (Delta_Eabs > tol) %  Check if Currents is non-zero
                K_I_P(l3,n3)  = YYL(i,j) .* (K_E.P(i,n3)  - K_E.P(j,n3))  + YYT(i,j) .* K_E.P(i,n3) ;
                K_I_Q(l3,n3)  = YYL(i,j) .* (K_E.Q(i,n3)  - K_E.Q(j,n3))  + YYT(i,j) .* K_E.Q(i,n3) ;
                K_I_Eabs(l3,n3) = YYL(i,j) .* (K_E.Eabs(i,n3) - K_E.Eabs(j,n3)) + YYT(i,j) .* K_E.Eabs(i,n3) ;
                K_I_Eang(l3,n3) = YYL(i,j) .* (K_E.Eang(i,n3) - K_E.Eang(j,n3)) + YYT(i,j) .* K_E.Eang(i,n3) ;

                K_Iabs_P(l3,n3) = (1 ./ abs(I(l3))) .* real(conj(I(l3)) .* K_I_P(l3,n3));
                K_Iabs_Q(l3,n3) = (1 ./ abs(I(l3))) .* real(conj(I(l3)) .* K_I_Q(l3,n3));
                K_Iabs_Eabs(l3,n3) = (1 ./ abs(I(l3))) .* real(conj(I(l3)) .* K_I_Eabs(l3,n3));
                K_Iabs_Eang(l3,n3) = (1 ./ abs(I(l3))) .* real(conj(I(l3)) .* K_I_Eang(l3,n3));

            else
                K_I_P(l3,n3) = zeros(nph, nph);
                K_I_Q(l3,n3) = zeros(nph, nph);
                K_I_Eabs(l3,n3) = zeros(nph, nph);
                K_I_Eang(l3,n3) = zeros(nph, nph);
                
                K_Iabs_P(l3,n3) = zeros(nph, nph);
                K_Iabs_Q(l3,n3) = zeros(nph, nph);
                K_Iabs_Eabs(l3,n3) = zeros(nph, nph);
                K_Iabs_Eang(l3,n3) = zeros(nph, nph);

            end
        end
    end
    
    %needed for the DC, otherwise if I -> 0, it messes up completely
%     K_Iabs_P(end-6:end-3,:) = -real(K_I_P(end-6:end-3,:)); 
%     K_Iabs_Q(end-6:end-3,:) = -real(K_I_Q(end-6:end-3,:)); 
%     K_Iabs_Eabs(end-6:end-3,:) = -real(K_I_Eabs(end-6:end-3,:)); 
%     K_Iabs_Eang(end-6:end-3,:) = -real(K_I_Eang(end-6:end-3,:)); 
%     
%     K_Iabs_P(end-2:end,:) = -real(K_I_P(end-2:end,:)); 
%     K_Iabs_Q(end-2:end,:) = -real(K_I_Q(end-2:end,:)); 
%     K_Iabs_Eabs(end-2:end,:) = -real(K_I_Eabs(end-2:end,:)); 
%     K_Iabs_Eang(end-2:end,:) = -real(K_I_Eang(end-2:end,:)); 
    
end
