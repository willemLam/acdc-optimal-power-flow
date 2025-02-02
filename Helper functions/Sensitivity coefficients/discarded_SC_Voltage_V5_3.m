function [K, Time] = SC_Voltage_V5_3(S_star,E_star,idx1,idx3,Grid_para,Filter_para,idxCtrl,unblanced_3ph,filter)

Yac = Grid_para.Yac;
Ydc = Grid_para.Ydc;
Sac = S_star(1:Grid_para.n_ph*Grid_para.n_ac);
Eac = E_star(1:Grid_para.n_ph*Grid_para.n_ac);
Sdc = S_star(Grid_para.n_ph*Grid_para.n_ac+1:end);
Edc = E_star(Grid_para.n_ph*Grid_para.n_ac+1:end);
n_ph = Grid_para.n_ph;
Fl = Grid_para.pos_ac3;

% the function computes analytically the SC of the nodes when one or more
% VSC's are present that operate on QV control. The function takes the filter into account

% This function computes the nodal voltage sensitivity coefficients for
% each control variable of the nodes listed in idxCtrl
% INPUT
% - Y           nodal admittance matrix
% - S0          apparent power injections at nominal voltage for all buses
% - idx1ph      structure containing network indices of all nodes  
%    .slack    indices of the slack buses
%    .pq       indices of the PQ buses
%    .qv       indices of the QV buses
% - idx3ph      structure containing 3ph expanded indices of all nodes
%    .slack    indices of the slack buses
%    .pq       indices of the PQ buses
%    .qv       indices of the QV buses
% - idxCtrl     vector containing netowrk indices of the nodes for which
%               the user wants to compute the voltage sensitivity
%               coefficients (SCs)
% - nph         number of phases in considered grid
% - vdep        structure containing the voltage-dependent weights and
%               exponents
%    .alpha    weights active power voltage-dependency
%    .lambda    exponents active power voltage-dependency
%    .beta     weights reactive power voltage-dependency
%    .omega    exponents reactive power voltage-dependency
% - tol         tolerance for Newton-Raphson convergence criterion
% - n_max       maximum number of iterations
%
% OUTPUT
% - K           A cell structure with complex, magnitude and angle voltage SCs
%               the size is length(idxCtrl) x 4, where
%   - column 1, has the index of the node
%   - column 2 a 2x1 cell containing each the complex voltage SCs
%   pertaining to the control variables of the node. 
%   . If a node is PQ then the 1st cell contains SCs w.r.t. P_{0,l}^{\phi} 
%   (active power) and the 2nd cell contains SCs w.r.t. Q_{0,l}^{\phi} (reactive power). 
%   . If a node is QV then the 1st cell contains SCs w.r.t. P_{0,m}^{\phi} (active power) and the
%   2nd cell contains SCs w.r.t. |\bar{E}^{\phi}_m| (voltage magnitude). 
%   . If a node is slack then the 1st cell contains SCs w.r.t. |\bar{E}^{\phi}_k| (voltage magnitude) and the
%   2nd cell contains SCs w.r.t. \angle(\bar{E}^{\phi}_k) (phase-angle).
%   (!!!) Each cell contains a nph|N| x nph matrix where |N| is the number of nodes (1ph) and
%   each column corresponds to the SCs pertaining to a phase and a control variable as explained above.
%   - column 3&4 have the same structure as column 2 but contain,
%   respectively, the magnitude and phase-angle nodal voltage SCs.


% Y = YY_augmented3;
% S0 = S_star;
% E = E_star;




% Bcell = mat2cell(B, 3*ones(1,15), 3*ones(1,15));
% A = arrayfun(@(x) A_inv*cell2mat(x), Bcell, 'UniformOutput', false)
% A = cell2mat(A)


%% Construct A
% Timing
T = tic;

if unblanced_3ph
    alp = exp(2*pi/3*1i);
    A_inv = (1/3*[1 1     1; 
         1 alp   alp^2; 
         1 alp^2 alp]);
    b_inv = [0;1;0];
else
    A_inv = eye(n_ph);
    b_inv = ones(n_ph,1);
end

if n_ph == 1
    idx = idx1;
else
    idx = idx3;
end



% X_mat = repmat({[ones(3,1),zeros(3,1),zeros(3,1)]}, 1, 2*(size(Eac,1)/nph)-2);
% X_mat = blkdiag(X_mat{:},eye(4));


A_invl = repmat({A_inv}, 1, size(Eac,1)/n_ph);
A_invlm = blkdiag(A_invl{:}); %cell2mat(A_invl); %

b_inv4 = repmat(b_inv,size(Fl,1)/n_ph,1);

for i = 1:n_ph:size(Yac,1)
    for j = 1:n_ph:size(Yac,2)
        Yac_s(i:i+n_ph-1,j:j+n_ph-1) = A_inv.*(transpose(sum(Yac(i:i+n_ph-1,j:j+n_ph-1))) );
    end
end
    
% Create F_{in}^{\phi\phi'}, P0 and Q0
Fac = diag(Eac)*conj(Yac)*diag(conj(Eac));
Fac_s = (A_invlm.*Eac).*conj(Yac)*(conj(A_invlm.*Eac));% diag(A_invlm*Eac)*conj(Yac)*diag(conj(A_invlm*Eac)); %diag(A_invlm*Eac)*conj(A_invlm*(Yac)*diag(Eac));%
Fac_s = diag(Eac)*conj(Yac_s)*diag(conj(Eac));

for i = 1:n_ph:size(Fac,1)
    for j = 1:n_ph:size(Fac,2)
        Fac_s(i:i+n_ph-1,j:j+n_ph-1) = A_inv.*(transpose(sum(Fac(i:i+n_ph-1,j:j+n_ph-1))) );
    end
end

P0ac = real(Sac);
Q0ac = imag(Sac);

% Fdc = diag(Edc)*conj(Ydc(19:26,19:26))*diag(conj(Edc));
Fdc = real(diag(Edc)*conj(Ydc)*diag(conj(Edc)));
P0dc = real(Sdc);
Q0dc = imag(Sdc);
      
F = blkdiag(Fac,Fdc);
F_s = blkdiag(Fac_s,Fdc);
E = [Eac;Edc];
E_s = [Eac;Edc];


%% Filter
I_b =Grid_para.Y_b*Grid_para.V_b;
Imag = abs(Grid_para.YY*E_star);
R_eq_ctu = interp1(Filter_para.IGBT_piecewise(:,1),Filter_para.IGBT_piecewise(:,2),Imag*I_b)./Grid_para.V_b./Imag;
R_eq = (R_eq_ctu*4/pi); % ???
R_eq(isnan(R_eq))=0;
R_eq(isinf(R_eq))=0;
Zf = (Filter_para.R + R_eq(sort([idx3.vscac_pq;idx3.vscac_vq]))) + 1i*Filter_para.X;

    
% influence of the filter

% alpha = ((conj(diag(Yac(Fl(:,1),Fl(:,2)))).*conj(Eac(Fl(:,2))) + conj(diag(Yac(Fl(:,1),Fl(:,1) ))).*conj(Eac(Fl(:,1)))).*conj(Zf(1:3:end)).*diag(Yac(Fl(:,1),Fl(:,1))).*Eac(Fl(:,1)) +  ...
%         ((diag(Yac(Fl(:,1),Fl(:,2)))).*(Eac(Fl(:,2))) + (diag(Yac(Fl(:,1),Fl(:,1) ))).*(Eac(Fl(:,1)))).*conj(Zf(1:3:end)).*conj(diag(Yac(Fl(:,1),Fl(:,1)))).*conj(Eac(Fl(:,1))) )./abs(Eac(Fl(:,1)));
% alpha = alpha(1:size(Fl,1));

% beta = ((conj(diag(Yac(Fl(:,1),Fl(:,2)))).*conj(Eac(Fl(:,2))) + conj(diag(Yac(Fl(:,1),Fl(:,1) ))).*conj(Eac(Fl(:,1)))).*conj(Zf(1:3:end)).*diag(Yac(Fl(:,1),Fl(:,1))).*Eac(Fl(:,1)) -  ...
%         ((diag(Yac(Fl(:,1),Fl(:,2)))).*(Eac(Fl(:,2))) + (diag(Yac(Fl(:,1),Fl(:,1) ))).*(Eac(Fl(:,1)))).*conj(Zf(1:3:end)).*conj(diag(Yac(Fl(:,1),Fl(:,1) ))).*conj(Eac(Fl(:,1))) );
% beta = beta(1:size(Fl,1));

if filter
    alpha = (conj(Yac(Fl(:,1),Fl(:,2)))*conj(Eac(Fl(:,2))) + conj(Yac(Fl(:,1),Fl(:,1)))*conj(Eac(Fl(:,1)))).*conj(Zf);
    beta = (Yac(Fl(:,1),Fl(:,2))*Eac(Fl(:,2)) + Yac(Fl(:,1),Fl(:,1))*Eac(Fl(:,1))).*conj(Zf);
    a = real(alpha.*(Yac(Fl(:,1),Fl(:,2))*Eac(Fl(:,2))) + beta.*(conj(Yac(Fl(:,1),Fl(:,2)))*conj(Eac(Fl(:,2)))));
    b = real(alpha.*(Yac(Fl(:,1),Fl(:,1))*Eac(Fl(:,1))) + beta.*(conj(Yac(Fl(:,1),Fl(:,1)))*conj(Eac(Fl(:,1)))));
    c = -imag(alpha.*(Yac(Fl(:,1),Fl(:,2))*Eac(Fl(:,2))) + beta.*(conj(Yac(Fl(:,1),Fl(:,2)))*conj(Eac(Fl(:,2)))));
    d = -imag(alpha.*(Yac(Fl(:,1),Fl(:,1))*Eac(Fl(:,1))) + beta.*(conj(Yac(Fl(:,1),Fl(:,1)))*conj(Eac(Fl(:,1)))));
    aa = zeros(length(idx.vscac_vq),length(idx.pqac)); aa(sub2ind(size(aa), (1:size(Fl,1))',Fl(:,2)-n_ph)) = a;
    bb = zeros(length(idx.vscac_vq),length(idx.vscac_vq)); bb(sub2ind(size(bb), (1:size(Fl,1))',(1:size(Fl,1))')) = b;
    cc = zeros(length(idx.vscac_vq),length(idx.pqac)); cc(sub2ind(size(cc), (1:size(Fl,1))',Fl(:,2)-n_ph)) = c;
    dd = zeros(length(idx.vscac_vq),length(idx.vscac_vq)); dd(sub2ind(size(dd), (1:size(Fl,1))',(1:size(Fl,1))')) = d;
else
    aa = 0;
    bb = 0;
    cc = 0;
    dd = 0;
end


% idx3ph.vscdc_vq = idx3ph.vscdc_vq-n_nodesac;
% idx3ph.pdc  =idx3ph.pdc-n_nodesac;
% Create each sub-matrix


%% PQac real
% A11
tmp_1 = real(F)./repmat(abs(E).',[size(F,1) 1]); tmp_1 = tmp_1(idx.pqac,idx.pqac);
tmp_3 = sum(real(F),2)./abs(E); tmp_3 = tmp_3(idx.pqac);
A11 = tmp_1 + diag(tmp_3);
% A12
tmp_1 = imag(F); tmp_1 = tmp_1(idx.pqac,idx.pqac);
tmp_2 = -sum(imag(F),2); tmp_2 = tmp_2(idx.pqac);
A12 = tmp_1 + diag(tmp_2);
% A13
tmp_1 = real(F)./repmat(abs(E).',[size(F,1) 1]); tmp_1 = tmp_1(idx.pqac,idx.vscac_vq);
A13 = tmp_1;

tmp_1 = imag(F); tmp_1 = tmp_1(idx.pqac,idx.vscac_vq);
A14 = tmp_1;

A15 = zeros(size(idx.pqac,1), size(idx.pdc,1));


%% PQac imag
% A21
tmp_1 = imag(F)./repmat(abs(E).',[size(F,1) 1]); tmp_1 = tmp_1(idx.pqac,idx.pqac);
tmp_3 = sum(imag(F),2)./abs(E); tmp_3 = tmp_3(idx.pqac);
A21 = tmp_1 + diag(tmp_3);
% A22
tmp_1 = -real(F); tmp_1 = tmp_1(idx.pqac,idx.pqac);
tmp_2 = sum(real(F),2); tmp_2 = tmp_2(idx.pqac);
A22 = tmp_1 + diag(tmp_2);
% A23
tmp_1 = imag(F)./repmat(abs(E).',[size(F,1) 1]); tmp_1 = tmp_1(idx.pqac,idx.vscac_vq);
A23 = tmp_1;

tmp_1 = -real(F); tmp_1 = tmp_1(idx.pqac,idx.vscac_vq);
A24 = tmp_1;

A25 = zeros(size(idx.pqac,1), size(idx.pdc,1));



%% VSC Pac = Pdc 
% A31
tmp_1 = real(F_s)./repmat(abs(E_s).',[size(F,1) 1]); tmp_1 = tmp_1(idx.vscac_vq,idx.pqac);
tmp_3 = sum(real(F_s),2)./abs(E_s); tmp_3 = tmp_3(idx.pqac);
% tmp_4 = zeros(length(idx3ph.pqac),1); tmp_4(Fl(:,1)) = real(alpha(1));
A31 = tmp_1 + aa;% + diag(tmp_3);
% A12
tmp_1 = imag(F_s); tmp_1 = tmp_1(idx.vscac_vq,idx.pqac);
tmp_2 = sum(imag(F_s),2); tmp_2 = tmp_2(idx.pqac);
A32 = tmp_1 + cc;% + diag(tmp_2)*0;
% A33
tmp_1 = real(F_s)./repmat(abs(E_s).',[size(F,1) 1]); tmp_1 = tmp_1(idx.vscac_vq,idx.vscac_vq);
tmp_3 = sum(real(F_s),2)./abs(E_s); tmp_3 = tmp_3(idx.vscac_vq);
A33 = tmp_1 + diag(tmp_3) + bb;
% A34
tmp_1 = imag(F_s); tmp_1 = tmp_1(idx.vscac_vq,idx.vscac_vq);
tmp_2 = -sum(imag(F_s),2); tmp_2 = tmp_2(idx.vscac_vq);
A34 = tmp_1 + diag(tmp_2) + dd;
% A35
% tmp_1 = real(F)./repmat(abs(E).',[size(F,1) 1]); tmp_1 = tmp_1(idx3ph.pdc,idx3ph.pdc);
% tmp_3 = sum(real(F),2)./abs(E); tmp_3 = tmp_3(idx3ph.vscdc_vq);
% A35 = repmat(-tmp_1 + -diag(tmp_3),[nph,1]); %Added - - 

% tmp_1 = real(F)./repmat(abs(E).',[size(F,1) 1]); tmp_1 = tmp_1(idx3ph.vscdc_vq,idx3ph.pdc);
% tmp_3 = sum(real(F),2)./abs(E); tmp_3 = tmp_3(idx3ph.vscdc_vq);
% A35 = -tmp_1;%-repmat(-tmp_1 + -diag(tmp_3),[nph,1]); %Added - - 
tmp_1 = real(F)./repmat(abs(E).',[size(F,1) 1]); tmp_1 = tmp_1(idx.vscdc_vq,idx.pdc); %% I changed this one
tmp_3 = sum(real(F),2)./abs(E); tmp_3 = tmp_3(idx.vscdc_vq);
tmp_4 = repmat(diag(tmp_1) + tmp_3,1,n_ph);
A35 = blkdiag(tmp_4(1,:)'.*b_inv,tmp_4(2,:)'.*b_inv,tmp_4(3,:)'.*b_inv,tmp_4(4,:)'.*b_inv);


%% VSC Q
% A41
tmp_1 = imag(F_s)./repmat(abs(E_s).',[size(F,1) 1]); tmp_1 = tmp_1(idx.vscac_vq,idx.pqac);
tmp_3 = sum(imag(F_s),2)./abs(E_s); tmp_3 = tmp_3(idx.pqac);
A41 = tmp_1;% + diag(tmp_3);
% A42
tmp_1 = -real(F_s); tmp_1 = tmp_1(idx.vscac_vq,idx.pqac);
tmp_2 = sum(real(F_s),2); tmp_2 = tmp_2(idx.pqac);
A42 = tmp_1;% + diag(tmp_2);
% A43
tmp_1 = imag(F_s)./repmat(abs(E_s).',[size(F,1) 1]); tmp_1 = tmp_1(idx.vscac_vq,idx.vscac_vq);
tmp_3 = sum(imag(F_s),2)./abs(E_s); tmp_3 = tmp_3(idx.vscac_vq);
A43 = tmp_1 + diag(tmp_3);
% A44
tmp_1 = -real(Fac_s); tmp_1 = tmp_1(idx.vscac_vq,idx.vscac_vq);
tmp_2 = sum(real(Fac_s),2); tmp_2 = tmp_2(idx.vscac_vq);
A44 = tmp_1 + diag(tmp_2);
% A45
A45 = zeros(size(idx.vscac_vq,1), size(idx.pdc,1));


%% PQdc
% A51
A51 = zeros(size(idx.pdc,1), size(idx.pqac,1));
A52 = zeros(size(idx.pdc,1), size(idx.pqac,1));
A53 = zeros(size(idx.pdc,1), size(idx.vscac_vq,1));
A54 = zeros(size(idx.pdc,1), size(idx.vscac_vq,1));
% A55
tmp_1 = real(F)./repmat(abs(E).',[size(F,1) 1]); tmp_1 = tmp_1(idx.pdc,idx.pdc);
tmp_3 = sum(real(F),2)./abs(E); tmp_3 = tmp_3(idx.pdc);
A55 = tmp_1 + diag(tmp_3); 



% A_inv = eye(3);
% if nph == 1
%  A_inv = 1;
% end
alp = exp(2*pi/3*1i);
    A_inv = 1/3*[1 1     1; 
             1 alp   alp^2; 
             1 alp^2 alp];
% if unblanced_3ph
%     for i = 1:3:size(A31,1)
%         for j = 1:3:size(A31,2)
%             A31(i:i+2,j:j+2) = real(A_inv.*(transpose(sum(A31(i:i+2,j:j+2))) + 1i*transpose(sum(A32(i:i+2,j:j+2)))));
%             A32(i:i+2,j:j+2) = imag(A_inv.*(transpose(sum(A31(i:i+2,j:j+2))) + 1i*transpose(sum(A32(i:i+2,j:j+2)))));
%         end
%     end
%     
%     
%     
% end
% if unblanced_3ph
%     
%     alp = exp(2*pi/3*1i);
%     A_inv = 1/3*[1 1     1; 
%              1 alp   alp^2; 
%              1 alp^2 alp];
%     A_inv4 = repmat({A_inv}, 1, size(Fl,1)/nph);
%     b_inv = [0;1;0];
%     b_inv4 = repmat(b_inv,size(Fl,1)/nph,1);
%     
%     A31((1:size(Fl,1))',Fl(:,2)-nph) = real(blkdiag(A_inv4{:})*complex(diag(A31((1:size(Fl,1))',Fl(:,2)-nph)),diag(A41((1:size(Fl,1))',Fl(:,2)-nph))).*blkdiag(A_inv4{:}));
%     A41((1:size(Fl,1))',Fl(:,2)-nph) = imag(blkdiag(A_inv4{:})*complex(diag(A31((1:size(Fl,1))',Fl(:,2)-nph)),diag(A41((1:size(Fl,1))',Fl(:,2)-nph))).*blkdiag(A_inv4{:}));
%     
%     A32((1:size(Fl,1))',Fl(:,2)-nph) = real(blkdiag(A_inv4{:})*complex(diag(A32((1:size(Fl,1))',Fl(:,2)-nph)),diag(A42((1:size(Fl,1))',Fl(:,2)-nph))).*blkdiag(A_inv4{:}));
%     A42((1:size(Fl,1))',Fl(:,2)-nph) = imag(blkdiag(A_inv4{:})*complex(diag(A32((1:size(Fl,1))',Fl(:,2)-nph)),diag(A42((1:size(Fl,1))',Fl(:,2)-nph))).*blkdiag(A_inv4{:}));
% 
%     A33((1:size(Fl,1))',(1:size(Fl,1))') = real(blkdiag(A_inv4{:})*complex(diag(A33((1:size(Fl,1))',(1:size(Fl,1))')),diag(A43((1:size(Fl,1))',(1:size(Fl,1))'))).*blkdiag(A_inv4{:}));
%     A43((1:size(Fl,1))',(1:size(Fl,1))') = imag(blkdiag(A_inv4{:})*complex(diag(A33((1:size(Fl,1))',(1:size(Fl,1))')),diag(A43((1:size(Fl,1))',(1:size(Fl,1))'))).*blkdiag(A_inv4{:}));
%     
%     A34((1:size(Fl,1))',(1:size(Fl,1))') = real(blkdiag(A_inv4{:})*complex(diag(A34((1:size(Fl,1))',(1:size(Fl,1))')),diag(A44((1:size(Fl,1))',(1:size(Fl,1))'))).*blkdiag(A_inv4{:}));
%     A44((1:size(Fl,1))',(1:size(Fl,1))') = imag(blkdiag(A_inv4{:})*complex(diag(A34((1:size(Fl,1))',(1:size(Fl,1))')),diag(A44((1:size(Fl,1))',(1:size(Fl,1))'))).*blkdiag(A_inv4{:}));
%     
% %     A31((1:size(Fl,1))',Fl(:,2)-nph) = real(blkdiag(A_inv4{:}).*complex(diag(A31((1:size(Fl,1))',Fl(:,2)-nph)),diag(A41((1:size(Fl,1))',Fl(:,2)-nph))));
% %     A41((1:size(Fl,1))',Fl(:,2)-nph) = imag(blkdiag(A_inv4{:}).*complex(diag(A31((1:size(Fl,1))',Fl(:,2)-nph)),diag(A41((1:size(Fl,1))',Fl(:,2)-nph))));
% %     
% %     A32((1:size(Fl,1))',Fl(:,2)-nph) = real(blkdiag(A_inv4{:}).*complex(diag(A32((1:size(Fl,1))',Fl(:,2)-nph)),diag(A42((1:size(Fl,1))',Fl(:,2)-nph))));
% %     A42((1:size(Fl,1))',Fl(:,2)-nph) = imag(blkdiag(A_inv4{:}).*complex(diag(A32((1:size(Fl,1))',Fl(:,2)-nph)),diag(A42((1:size(Fl,1))',Fl(:,2)-nph))));
% % 
% %     A33((1:size(Fl,1))',(1:size(Fl,1))') = real(blkdiag(A_inv4{:}).*complex(diag(A33((1:size(Fl,1))',(1:size(Fl,1))')),diag(A43((1:size(Fl,1))',(1:size(Fl,1))'))));
% %     A43((1:size(Fl,1))',(1:size(Fl,1))') = imag(blkdiag(A_inv4{:}).*complex(diag(A33((1:size(Fl,1))',(1:size(Fl,1))')),diag(A43((1:size(Fl,1))',(1:size(Fl,1))'))));
% %     
% %     A34((1:size(Fl,1))',(1:size(Fl,1))') = real(blkdiag(A_inv4{:}).*complex(diag(A34((1:size(Fl,1))',(1:size(Fl,1))')),diag(A44((1:size(Fl,1))',(1:size(Fl,1))'))));
% %     A44((1:size(Fl,1))',(1:size(Fl,1))') = imag(blkdiag(A_inv4{:}).*complex(diag(A34((1:size(Fl,1))',(1:size(Fl,1))')),diag(A44((1:size(Fl,1))',(1:size(Fl,1))'))));
% 
%     
%     tmp_4 = reshape(b_inv4 .* A35(A35~=0),nph,size(Fl,1)/nph)';
%     A35 = blkdiag(tmp_4(1,:)',tmp_4(2,:)',tmp_4(3,:)',tmp_4(4,:)');
% 
% end
%% Assemble A
A = [ A11 A12 A13 A14 A15; ...
      A21 A22 A23 A24 A25; ...
      A31 A32 A33 A34 A35; ...
      A41 A42 A43 A44 A45; ...
      A51 A52 A53 A54 A55];

Time.A = toc(T);

%%
% 
% %% eq 3
% 
% sc = complex(J_VR(:,1),J_VX(:,1))
% sc_mag = 1./abs(E(2:end)) .* real(conj(E(2:end)) .* sc)
% sc_ang = 1./abs(E(2:end)).^2 .* imag(conj(E(2:end)) .* sc)
% 
% sum(sum(real(F(15,:)))/abs(E(15)) + real(F(15,15))/abs(E(15)))*sc_mag(15-1)+...  %15
% sum(real(F(15,[1:14,16]))./transpose(abs(E([1:14,16]))))*sc_mag(7-1)+...  %1:14,16
% sum(-sum(imag(F(15,:))) - imag(F(15,15)))*sc_ang(15-1)+...  %15
% sum(imag(F(15,[1:14,16])))*sc_ang(7-1)+... %1:14,16
% sum(-real(F(17,[18:20]))./transpose(abs(E([18:20]))))*sc_mag(19-1) %18:20
% 
% sum(sum(real(F(17,:)))/abs(E(17)) + real(F(17,17))/abs(E(17)))*sc_mag(17-1) %17
% 
% 
% %validated -> eq3 is correct
% 
% sum(sum(real(F(15,:)))/abs(E(15)) + real(F(15,15))/abs(E(15)))*J_VR(15-1,1) +... %15
% sum(real(F(15,[1:14,16]))./transpose(abs(E([1:14,16]))))*J_VR(7-1,1) +...%1:14,16
% sum(-sum(imag(F(15,:))) - imag(F(15,15)))*J_VX(15-1,1) +...%15
% sum(imag(F(15,[1:14,16])))*J_VX(7-1,1) %1:14,16
% 
% sum(sum(real(F(17,:)))/abs(E(17)) + real(F(17,17))/abs(E(17)))*J_VR(17-1,1) +...%17
% sum(-real(F(17,[18:20]))./transpose(abs(E([18:20]))))*J_VR(19-1,1) %18:20
% 
% 
% sum(real(F(16,:)))/abs(E(16)) + real(F(16,16))/abs(E(16)) %16
% real(F(16,[1:14,15]))./transpose(abs(E([1:14,15]))) %1:14,15
% 
% sum(real(F(18,:)))/abs(E(18)) + real(F(18,18))/abs(E(18)) %18
% -real(F(18,[17,19:20]))./transpose(abs(E([17,19:20]))) %17,19:20
% 
% -sum(imag(F(16,:))) - imag(F(16,16)) %16
% imag(F(16,[1:14,15])) %1:14,15
% 
% %% eq 4
% 
% sc = complex(J_QR(:,14),J_QX(:,14))
% sc_mag = 1./abs(E(2:end)) .* real(conj(E(2:end)) .* sc)
% sc_ang = 1./abs(E(2:end)).^2 .* imag(conj(E(2:end)) .* sc)
% 
% (sum(imag(F(15,:)))/abs(E(15)) + imag(F(15,15))/abs(E(15)))*sc_mag(15-1)+... %15
% sum(imag(F(15,[1:14,16]))./transpose(abs(E([1:14,16]))))*sc_mag(7-1)+... %1:14,16
% (sum(real(F(15,:))) - real(F(15,15)))*sc_ang(15-1)+... %15
% sum(-real(F(15,[1:14,16])))*sc_ang(7-1) %1:14,16
% 
% %validated -> sum = -1
% A4 = [ A41 A42 A43 A44 A45]
% %%
% 
% %% eq 5
% 
% sc = complex(J_PR(:,14),J_PX(:,14))
% sc_mag = 1./abs(E(2:end)) .* real(conj(E(2:end)) .* sc)
% sc_ang = 1./abs(E(2:end)).^2 .* imag(conj(E(2:end)) .* sc)
% 
% (sum(real(F(19,:)))/abs(E(19)) + real(F(19,19))/abs(E(19)))*sc_mag(19-1)+... %15
% sum(real(F(19,[17,18,20]))./transpose(abs(E([17,18,20])))).*sc_mag([17,18,20]-1) %1:14,16
% 
% %validated -> sum = 1
% A5 = [ A51 A52 A53 A54 A55]



 %% Compute nodal voltage SCs for each control variable

% Initialize Outputs
K = cell(length(idxCtrl),4);
Time.K = [];

for id_x = 1:length(idxCtrl)
    % Initialize cell entry
    K{id_x,1} = idxCtrl(id_x); % 1ph index
    K{id_x,2} = cell(2,1);     % Nodal Voltage SCs (complex)
    K{id_x,3} = cell(2,1);     % Magnitude Voltage SCs
    K{id_x,4} = cell(2,1);     % Angle Voltage SCs
    
    % Identify node-type of idxCtrl [1: slack, 2: pq, 3: qv]
     node_type = 1*sum( idxCtrl(id_x) == idx.slack ) + ...
                 2*sum( idxCtrl(id_x) == idx.pqac ) + ...
                 3*sum( idxCtrl(id_x) == idx.vscac_vq ) + ...
                 4*sum( idxCtrl(id_x) == idx.vscdc_vq ) + ...
                 5*sum( idxCtrl(id_x) == idx.pdc ); 
    for ctrl_var = 1:2
       
     T = tic;   
           
     % Initialize outputs
     K{id_x,2}{ctrl_var,1} = zeros(size(E,1),1); % complex
     K{id_x,3}{ctrl_var,1} = zeros(size(E,1),1); % magnitude
     K{id_x,4}{ctrl_var,1} = zeros(size(E,1),1); % angle

     % Get 3ph index
     tmp_idx = idxCtrl(id_x);
     
     % Compute magnitude and angle nodal voltage SCs for slack nodes
     switch( node_type )
         case 1 % Slack node
             if(ctrl_var == 1) % if X = |\bar{E}^{\phi}_k|
                %K{id_x,2}{ctrl_var,1}(tmp_idx,nph) = exp(1i*angle(E(tmp_idx,1))); % complex
                K{id_x,3}{ctrl_var,1}(tmp_idx,1) = 1; % magnitude
             else % if X = \angle(\bar{E}^{\phi}_k)
                %K{id_x,2}{ctrl_var,1}(tmp_idx,nph) = 1i*E(tmp_idx,1); % complex
                K{id_x,4}{ctrl_var,1}(tmp_idx,1) = 1; % angle
             end
         otherwise % non-slack node
             % Do nothing as derivatives are all zeros.
     end
     
     % Compute magnitude nodal voltage SCs for QV nodes
     switch( node_type )
         case 4 % vscdc_vq node
             if(ctrl_var == 1) % if X = |\bar{E}^{\phi}_m|
                K{id_x,3}{ctrl_var,1}(tmp_idx,1) = 1; % magnitude
             else % if X = P0^{\phi}_{0,m}
                % Do nothing as derivative is zero
             end
         otherwise % non-qv node
             % Do nothing as derivatives are all zeros.
     end
     
     % Compute nodal voltage SCs for PQ & QV nodes
     % Construct u(X)
      switch( node_type )
         case 1 % slack node
             if(ctrl_var == 1) % if X = |\bar{E}^{\phi}_k|
                u1 = -((real(Fac(idx.pqac,tmp_idx))/(abs(Eac(tmp_idx))^2))*real(conj(Eac(tmp_idx))*exp(1i*angle(Eac(tmp_idx)))) + ...
                     (imag(Fac(idx.pqac,tmp_idx))/(abs(Eac(tmp_idx))^2))*imag(conj(Eac(tmp_idx))*exp(1i*angle(Eac(tmp_idx)))) );
                u2 = (real(Fac(idx.pqac,tmp_idx))/(abs(Eac(tmp_idx))^2))*imag(conj(Eac(tmp_idx))*exp(1i*angle(Eac(tmp_idx)))) - ...
                     (imag(Fac(idx.pqac,tmp_idx))/(abs(Eac(tmp_idx))^2))*real(conj(Eac(tmp_idx))*exp(1i*angle(Eac(tmp_idx))));
                u3 = zeros(length(idx.vscac_vq),1);
                u4 = zeros(length(idx.vscac_vq),1);
                u5 = zeros(length(idx.pdc),1);
             else % if X = \angle(\bar{E}^{\phi}_k)
                u1 = -((real(Fac(idx.pqac,tmp_idx))/(abs(Eac(tmp_idx))^2))*real(conj(Eac(tmp_idx))*1i*Eac(tmp_idx)) + ...
                     (imag(Fac(idx.pqac,tmp_idx))/(abs(Eac(tmp_idx))^2))*imag(conj(Eac(tmp_idx))*1i*Eac(tmp_idx)));
                u2 = (real(Fac(idx.pqac,tmp_idx))/(abs(Eac(tmp_idx))^2))*imag(conj(Eac(tmp_idx))*1i*Eac(tmp_idx)) - ...
                     (imag(Fac(idx.pqac,tmp_idx))/(abs(Eac(tmp_idx))^2))*real(conj(Eac(tmp_idx))*1i*Eac(tmp_idx));
                u3 = zeros(length(idx.vscac_vq),1);
                u4 = zeros(length(idx.vscac_vq),1);
                u5 = zeros(length(idx.pdc),1);
             end
          case 2 % pqac node
             if(ctrl_var == 1) % if X = P0^{\phi}_{0,l}
                u1 = zeros(length(idx.pqac),1);
                u1(tmp_idx == idx.pqac,1) = 1; %sum((vdep.alpha(tmp_idx,:)),2);
                u2 = zeros(length(idx.pqac),1);  
                u3 = zeros(length(idx.vscac_vq),1);
                u4 = zeros(length(idx.vscac_vq),1);
                u5 = zeros(length(idx.pdc),1);        
             else % if X = Q0^{\phi}_{0,l}
                u1 = zeros(length(idx.pqac),1);
                u2 = zeros(length(idx.pqac),1);
                u2(tmp_idx == idx.pqac,1) = 1; %sum((vdep.beta(tmp_idx,:)),2);                
                u3 = zeros(length(idx.vscac_vq),1);
                u4 = zeros(length(idx.vscac_vq),1);
                u5 = zeros(length(idx.pdc),1);     
             end   
          case 3 % vscac_vq node
             if(ctrl_var == 1) % if X = P0^{\phi}_{0,m}
                u1 = zeros(length(idx.pqac),1);
                u2 = zeros(length(idx.pqac),1);
                u3 = zeros(length(idx.vscac_vq),1);
                u4 = zeros(length(idx.vscac_vq),1);
                u4(tmp_idx == idx.vscac_vq,1) = -1;
                u5 = zeros(length(idx.pdc),1);  
                
             else % if X = |\bar{E}^{\phi}_m|
                u1 = zeros(length(idx.pqac),1);
                u2 = zeros(length(idx.pqac),1);
                u3 = zeros(length(idx.vscac_vq),1);
                u4 = zeros(length(idx.vscac_vq),1);
                u5 = zeros(length(idx.pdc),1);         
             end   
             
          case 4 % vscdc_vq node
             if(ctrl_var == 1) % if X = P0^{\phi}_{0,m}
              
                u1 = zeros(length(idx.pqac),1);
                u2 = zeros(length(idx.pqac),1);
                 
                u3 = zeros(length(idx.vscac_vq),1);
                tmp_1 = real(F_s)./repmat(abs(E_s).',[size(F,1) 1]); tmp_1 = tmp_1(idx.vscdc_vq,idx.vscdc_vq);
                tmp_3 = sum(real(F_s),2)./abs(E_s); tmp_3 = tmp_3(idx.vscdc_vq);
                tmp_4 = repmat(-diag(tmp_1) + -(tmp_3),[1,n_ph])';
                u3(reshape(repmat(tmp_idx == idx.vscdc_vq,1,n_ph)',size(Fl,1),1),1) = tmp_4(:,tmp_idx - idx.vscdc_vq(1) +1);%.*[0;1;0];%tmp_idx-idx3ph.vscdc_vq(1)+1);                
                
                u4 = zeros(length(idx.vscac_vq),1);
                
                u5 = zeros(length(idx.vscdc_vq),1);
                tmp_1 = real(F)./repmat(abs(E).',[size(F,1) 1]); tmp_1 = tmp_1(idx.pdc,idx.vscdc_vq);
                tmp_3 = sum(real(F),2)./abs(E); tmp_3 = tmp_3(idx.pdc);
                tmp_4 = -diag(tmp_1) + -tmp_3;
                u5(tmp_idx == idx.vscdc_vq,1) = tmp_4(tmp_idx - idx.vscdc_vq(1) +1);%tmp_idx-idx3ph.vscdc_vq(1)+1);
         
             else % if X = |\bar{E}^{\phi}_m|
                u1 = zeros(length(idx.pqac),1);
                u2 = zeros(length(idx.pqac),1);
                u3 = zeros(length(idx.vscac_vq),1);
                u4 = zeros(length(idx.vscac_vq),1);
                u5 = zeros(length(idx.pdc),1);         
             end   

          case 5 % pdc node
             if(ctrl_var == 1) % if X = P0^{\phi}_{0,l}
                u1 = zeros(length(idx.pqac),1);
                u2 = zeros(length(idx.pqac),1);
                u3 = zeros(length(idx.vscac_vq),1);
                u4 = zeros(length(idx.vscac_vq),1);
                u5 = zeros(length(idx.pdc),1);  
                u5(tmp_idx == idx.pdc,1) = 1;  %Added -
             else % if X = Q0^{\phi}_{0,l}
                u1 = zeros(length(idx.pqac),1);
                u2 = zeros(length(idx.pqac),1);
                u3 = zeros(length(idx.vscac_vq),1);
                u4 = zeros(length(idx.vscac_vq),1);
                u5 = zeros(length(idx.pdc),1);     
              end   
          otherwise % "error"
              % Do nothing
      end    
     
     if unblanced_3ph
         u3 = b_inv4.*u3; %kron(eye(size(Fl,1)/nph),A_inv)*u3;
         u4 = reshape(sum(reshape(u4,3,4)).*b_inv,12,1); %b_inv4.*u4; %kron(eye(size(Fl,1)/nph),A_inv)*u4;
     end
     
     % Solve A*x(X)=u(X) 
     u = [u1;u2;u3;u4;u5];
     x = linsolve(A,u);
     
%      x_unbalanced = x;
%      u_unbalanced = u;
%      A_unbalanced = A;

%      x_balanced = x;
%      u_balanced = u;
%      A_balanced = A;
     
     % Assemble magnitude nodal voltage SCs
     K{id_x,3}{ctrl_var,1}(idx.pqac,1) = ...
         x( 1:length(idx.pqac) ); % PQac nodes
     K{id_x,3}{ctrl_var,1}(idx.vscac_vq,1) = ...
         x( 2*length(idx.pqac) +1: 2*length(idx.pqac) +length(idx.vscac_vq)  ); % vsc nodes
     K{id_x,3}{ctrl_var,1}(idx.pdc,1) = ...  %(idx3ph.pdc+18,ph) = ...
         x( 2*length(idx.pqac) + 2*length(idx.vscac_vq) +1 : 2*length(idx.pqac) + 2*length(idx.vscac_vq) +length(idx.pdc) ); % pdc nodes
     
     % Assemble phase-angle nodal voltage SCs
     K{id_x,4}{ctrl_var,1}(idx.pqac,1) = ...
         x( length(idx.pqac)+1 : 2*length(idx.pqac) ); % PQac nodes     
     K{id_x,4}{ctrl_var,1}(idx.vscac_vq,1) = ...
         x( 2*length(idx.pqac) +length(idx.vscac_vq) +1: 2*length(idx.pqac) + 2*length(idx.vscac_vq)   ); % VSCac nodes  
     
     % Infer nodal voltage SCs 
     K{id_x,2}{ctrl_var,1}(:,1) = ...
        E.*( (1./abs(E)).*K{id_x,3}{ctrl_var,1}(:,1) + ...
            1i*K{id_x,4}{ctrl_var,1}(:,1) );
     
     % Time
     Time.K = [Time.K; toc(T)];
     
    end
 end

end