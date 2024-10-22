function [E, S, IB_big,  IDXall,  COEFF, COEFFq, COEFFp_I_mat, COEFFq_I_mat, C_rp , C_rq, C_xp, C_xq] = ...
    LFandCoeff(Y,YYL, YYT,n_lines, n_nodes, H, pload, qload, pgen_opt, qgen_opt, Ab, Vb, current_idx)
v_slack = 1;
for ts = 1: H
[E(:,ts), S(:,ts)]= ZNR_Load_Flow(Y, v_slack, [0;pload(:,ts)], [0;qload(:,ts)], pgen_opt(:,ts), qgen_opt(:,ts), Ab, Vb );

nph =1;
   for i = 1:nph:size(YYL,1)
            
            for j = 1:nph:size(YYL,2)
                I_YYT(ceil(i/nph),ceil(j/nph)) = YYL(i:i+(nph-1),j:j+(nph-1))*(E(i:i+(nph-1),ts) - E(j:j+(nph-1),ts)) + ...
                                            YYT(i:i+(nph-1),j:j+(nph-1))*E(i:i+(nph-1),ts);
            end
   end
   
IDXall = find(I_YYT>0);
% I_LF(:,ts) = I_YYT(IDXall);

slack2 =1;
slack = 1;
alphap = 0;
alphaq = 0;
E0 =1;
% [COEFF1{ts},COEFFcomplex1,COEFFq1{ts},COEFFqcomplex1]=Coeffs_Voltage_Alpha(Y,S(:,ts),abs(E(:,ts)),slack2,nph,alphap,alphaq,E0);

[COEFF{ts},COEFFcomplex,COEFFq{ts},COEFFqcomplex]=Coeffs_Voltage_Alpha(Y,S(:,ts),E(:,ts),slack2,nph,alphap,alphaq,E0);

%ReConstruct Complex Voltage Coeffs
COEFFcomplex = COEFFcomplex(1:0.5*size(COEFFcomplex,1),:) ...
              + 1i*COEFFcomplex(0.5*size(COEFFcomplex,1)+1:end, :);
COEFFqcomplex = COEFFqcomplex(1:0.5*size(COEFFqcomplex,1),:) ...
              + 1i*COEFFqcomplex(0.5*size(COEFFqcomplex,1)+1:end, :); 

K_p{ts} = [COEFFcomplex(1:(nph*(slack-1) + 1 - 1),:); zeros(nph,size(COEFFcomplex,2)) ; COEFFcomplex(nph*(slack-1)+1:end,:)];
K_q{ts} = [COEFFqcomplex(1:(nph*(slack-1) + 1 - 1),:); zeros(nph,size(COEFFqcomplex,2)) ; COEFFqcomplex(nph*(slack-1)+1:end,:)];

[C_rp{ts} , C_rq{ts}, C_xp{ts}, C_xq{ts}] = Coeffs_Losses(Y, E(:,ts), K_p{ts}, K_q{ts});

% Lossp = diffP(ts-1,:)*C_rp(ts-1,:)' + diffQ(ts-1,:)*C_rq(ts-1,:)';
% Lossq = diffP(ts-1,:)*C_xp(ts-1,:)' + diffQ(ts-1,:)*C_xq(ts-1,:)';

% delLosspq(ts-1) = complex(Lossp(ts-1), Lossq(ts-1));


%         K_p{ts-1,1} = [zeros(nph,size(COEFF,2)); COEFF];
%         K_q{ts-1,1} = [zeros(nph,size(COEFF,2)); COEFF];        


% [Icoeff, Icoeff_complex, Icoeffq,Icoeffq_complex, Currents]=Coeffs_Currents_org(YYL,YYT,E(:,ts),K_p{ts},K_q{ts},nph);
[Icoeff, Icoeff_complex, Icoeffq,Icoeffq_complex, COEFFp_I_mat1, COEFFq_I_mat1, Currents] = Coeffs_Currents_orgIII(YYL,YYT,E(:,ts),K_p{ts},K_q{ts},nph, current_idx);

%   for i = 1:n_lines
%             
%             COEFFp_I_mat{ts}(:,i) = reshape(cell2mat(Icoeff{i}),n_nodes*n_nodes,1);
%             COEFFq_I_mat{ts}(:,i) = reshape(cell2mat(Icoeffq{i}),n_nodes*n_nodes,1);
%             
%   end 
IB_big(:,ts) = abs(I_YYT(current_idx));

COEFFp_I_mat{ts} = COEFFp_I_mat1;
COEFFq_I_mat{ts} = COEFFq_I_mat1;

  
end

end


















