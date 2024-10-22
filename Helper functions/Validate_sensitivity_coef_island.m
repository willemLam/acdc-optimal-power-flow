% %% main script for linear power system state estimation
% 
% clear all;
% close all;
% clc;
% 
% addpath(genpath(pwd))
% 
% % get para
% [Grid_para, Filter_para, idx1, idx3, constraints] = initialize_RT_island();
% idx = idx1; %single phase equivalent
% 
% % get data
% port_SE = 35005;
% [E_star, S_star, I_abs] = get_states_from_SE(port_SE, Grid_para);
% E_star(abs(E_star)<0.5) = 0.9375;
% 
% % balanced
% modes = {'P9';'Q9';'Q22';'E26';'P30';'E19_abs';'E19_angle';'E27'};
% % modes = {'E26'};

%% main script for linear power system state estimation

clear all;
close all;
clc;

addpath(genpath(pwd))


%% Get the EMTP measurements

[Grid_para, Filter_para, idx1, idx3, constraints] = initialize_RT_island();
% [Grid_para, Filter_para, idx1, idx3, constraints] = initialize_RT();
idx = idx1; %single phase equivalent
tol = 1e-6;
n_max = 20;



[E_star, S_star, I_star] = get_fake_states(idx, Grid_para,Filter_para);
S_Loss = sum(S_star);




data_all_Vabs = {};
data_real_Vabs = {};
data_imag_Vabs = {};

data_all_Vang = {};
data_real_Vang = {};
data_imag_Vang = {};

data_all_I = {};
data_real_I = {};
data_imag_I = {};

data_all_Loss = {};
data_real_Loss = {};
data_imag_Loss = {};

%% GO
modes = {'P9';'Q9';'Q22';'E26';'P30';'E19_abs';'E19_angle';'E27'};
for m=1:length(modes)
    
    mode = char(modes(m));
    
    E_star2 = E_star;
    S_star2 = S_star;
    
    ad  = 1e-6;
    switch mode
        case 'P30'
            S_star2(idx.pdc(3)) = S_star2(idx.pdc(3)) + ad;
        case 'E26'
            E_star2(idx.vscdc_vq(2)) = E_star2(idx.vscdc_vq(2)) + ad;
        case 'Q22'
            S_star2(idx.vscac_vq(2)) = S_star2(idx.vscac_vq(2)) + ad*1i;
        case 'Q9'
            S_star2(idx.pqac(9)) = S_star2(idx.pqac(9)) + ad*1i;
        case 'P9'
            S_star2(idx.pqac(9)) = S_star2(idx.pqac(9)) + ad;
        case 'E19_abs'
            E_star2(idx.vscac_vv(1)) = E_star2(idx.vscac_vv(1)) * (1-ad);
        case 'E19_angle'
            E_star2(idx.vscac_vv(1)) = abs(E_star2(idx.vscac_vv(1)))*(cos(angle(E_star2(idx.vscac_vv(1)))+ad) + 1i*sin(angle(E_star2(idx.vscac_vv(1)))+ad) );
        case 'E27'
            E_star2(idx.vdc(1)) = E_star2(idx.vdc(1)) + ad;
        
        otherwise
            disp('NOPE') 
    end
    
    
    if Grid_para.n_ph == 3
        warning('not implemented yet')
    elseif Grid_para.n_ph == 1
        [E_star2,J2,~] = NR_rectangularACDC_1ph_general_V2(Grid_para,Filter_para,S_star2,E_star2,Grid_para.E_0,idx,tol,n_max);
    end
    S_star2 = E_star2.*conj(Grid_para.YY*E_star2);
    I_star2 = get_Current_flow(E_star2,Grid_para);
    S_Loss2 = sum(S_star2);

%% Compute SC
idxCtrl = 1:Grid_para.n_nodes; %estimate all voltage sensitivity coefficients
[SC, ~] = SC_voltage_rectangular_V4(E_star,idx,Grid_para,idxCtrl);
[K.Eabs.P, K.Eabs.Q, K.Eabs.Eabs,  K.Eabs.Eang, K.Eang.P, K.Eang.Q, K.Eang.Eabs, K.Eang.Eang] = transform_K_polar(SC, Grid_para, idx);
[K.E.P, K.E.Q, K.E.Eabs, K.E.Eang] = transform_K_complex(SC, Grid_para, idx);
[K.Iabs.P, K.I.P , K.Iabs.Q, K.I.Q , K.Iabs.Eabs, K.I.Eabs , K.Iabs.Eang, K.I.Eang  ] = Coeffs_Currents_GPT(I_star,E_star,K.E,Grid_para) ;
[K.Rloss.P, K.Rloss.Q, K.Rloss.Eabs, K.Rloss.Eang, K.Xloss.P, K.Xloss.Q,  K.Xloss.Eabs, K.Xloss.Eang ] = Coeffs_Losses(E_star, K.E,Grid_para);


%% Analyse results
disp(mode)
switch mode
case 'P30'
    i = 30; 
    
    r = K.Eabs.P(:,i);
    lf = (abs(E_star)-abs(E_star2))/real(S_star(i)-S_star2(i));
    data_real_Vabs.num(m,:) = [sqrt(mean((real(r) - real(lf)).^2)), max(abs(max((real(r) - real(lf)))), abs(min((real(r) - real(lf))))) ];
    data_imag_Vabs.num(m,:) = [sqrt(mean((imag(r) - imag(lf)).^2)), max(abs(max((imag(r) - imag(lf)))), abs(min((imag(r) - imag(lf))))) ];
    data_all_Vabs.num(m,:) = r - lf;
    
    r = K.Eang.P(:,i);
    lf = (angle(E_star)-angle(E_star2))/real(S_star(i)-S_star2(i));
    data_real_Vang.num(m,:) = [sqrt(mean((real(r) - real(lf)).^2)), max(abs(max((real(r) - real(lf)))), abs(min((real(r) - real(lf))))) ];
    data_imag_Vang.num(m,:) = [sqrt(mean((imag(r) - imag(lf)).^2)), max(abs(max((imag(r) - imag(lf)))), abs(min((imag(r) - imag(lf))))) ];
    data_all_Vang.num(m,:) = r - lf;
    
    r = K.Iabs.P(:,i);
    lf = (abs(I_star)-abs(I_star2))/real(S_star(i)-S_star2(i));
    data_real_I.num(m,:) = [sqrt(mean((real(r) - real(lf)).^2)), max(abs(max((real(r) - real(lf)))), abs(min((real(r) - real(lf))))) ];
    data_imag_I.num(m,:) = [sqrt(mean((imag(r) - imag(lf)).^2)), max(abs(max((imag(r) - imag(lf)))), abs(min((imag(r) - imag(lf))))) ];
    data_all_I.num(m,:) = r - lf;

    r = complex(K.Rloss.P(i),K.Xloss.P(i));
    lf = (S_Loss-S_Loss2)/real(S_star(i)-S_star2(i));
    data_real_Loss.num(m,:) = [sqrt(mean((real(r) - real(lf)).^2)), max(abs(max((real(r) - real(lf)))), abs(min((real(r) - real(lf))))) ];
    data_imag_Loss.num(m,:) = [sqrt(mean((imag(r) - imag(lf)).^2)), max(abs(max((imag(r) - imag(lf)))), abs(min((imag(r) - imag(lf))))) ];
    data_all_Loss.num(m,:) = r - lf;
    
case 'E26'
    i = 26; 
    r = K.Eabs.Eabs(:,i);
    lf = (abs(E_star)-abs(E_star2))/real(E_star(i)-E_star2(i));
    data_real_Vabs.num(m,:) = [sqrt(mean((real(r) - real(lf)).^2)), max(abs(max((real(r) - real(lf)))), abs(min((real(r) - real(lf))))) ];
    data_imag_Vabs.num(m,:) = [sqrt(mean((imag(r) - imag(lf)).^2)), max(abs(max((imag(r) - imag(lf)))), abs(min((imag(r) - imag(lf))))) ];
    data_all_Vabs.num(m,:) = r - lf;
    
    r = K.Eang.Eabs(:,i);
    lf = (angle(E_star)-angle(E_star2))/real(E_star(i)-E_star2(i));
    data_real_Vang.num(m,:) = [sqrt(mean((real(r) - real(lf)).^2)), max(abs(max((real(r) - real(lf)))), abs(min((real(r) - real(lf))))) ];
    data_imag_Vang.num(m,:) = [sqrt(mean((imag(r) - imag(lf)).^2)), max(abs(max((imag(r) - imag(lf)))), abs(min((imag(r) - imag(lf))))) ];
    data_all_Vang.num(m,:) = r - lf;
    
    r = K.Iabs.Eabs(:,i);
    lf = (abs(I_star)-abs(I_star2))/real(E_star(i)-E_star2(i));
    data_real_I.num(m,:) = [sqrt(mean((real(r) - real(lf)).^2)), max(abs(max((real(r) - real(lf)))), abs(min((real(r) - real(lf))))) ];
    data_imag_I.num(m,:) = [sqrt(mean((imag(r) - imag(lf)).^2)), max(abs(max((imag(r) - imag(lf)))), abs(min((imag(r) - imag(lf))))) ];
    data_all_I.num(m,:) = r - lf;

    r = complex(K.Rloss.Eabs(i),K.Xloss.Eabs(i));
    lf = (S_Loss-S_Loss2)/real(E_star(i)-E_star2(i));
    data_real_Loss.num(m,:) = [sqrt(mean((real(r) - real(lf)).^2)), max(abs(max((real(r) - real(lf)))), abs(min((real(r) - real(lf))))) ];
    data_imag_Loss.num(m,:) = [sqrt(mean((imag(r) - imag(lf)).^2)), max(abs(max((imag(r) - imag(lf)))), abs(min((imag(r) - imag(lf))))) ];
    data_all_Loss.num(m,:) = r - lf;
        
    
    
case 'Q22'
    i = 22;
    
    r = K.Eabs.Q(:,i);
    lf = (abs(E_star)-abs(E_star2))/imag(S_star(i(1))-S_star2(i(1)));
    data_real_Vabs.num(m,:) = [sqrt(mean((real(r) - real(lf)).^2)), max(abs(max((real(r) - real(lf)))), abs(min((real(r) - real(lf))))) ];
    data_imag_Vabs.num(m,:) = [sqrt(mean((imag(r) - imag(lf)).^2)), max(abs(max((imag(r) - imag(lf)))), abs(min((imag(r) - imag(lf))))) ];
    data_all_Vabs.num(m,:) = r - lf;
    
    r = K.Eang.Q(:,i);
    lf = (angle(E_star)-angle(E_star2))/imag(S_star(i(1))-S_star2(i(1)));
    data_real_Vang.num(m,:) = [sqrt(mean((real(r) - real(lf)).^2)), max(abs(max((real(r) - real(lf)))), abs(min((real(r) - real(lf))))) ];
    data_imag_Vang.num(m,:) = [sqrt(mean((imag(r) - imag(lf)).^2)), max(abs(max((imag(r) - imag(lf)))), abs(min((imag(r) - imag(lf))))) ];
    data_all_Vang.num(m,:) = r - lf;
    
    r = K.Iabs.Q(:,i);
    lf = (abs(I_star)-abs(I_star2))/imag(S_star(i(1))-S_star2(i(1)));
    data_real_I.num(m,:) = [sqrt(mean((real(r) - real(lf)).^2)), max(abs(max((real(r) - real(lf)))), abs(min((real(r) - real(lf))))) ];
    data_imag_I.num(m,:) = [sqrt(mean((imag(r) - imag(lf)).^2)), max(abs(max((imag(r) - imag(lf)))), abs(min((imag(r) - imag(lf))))) ];
    data_all_I.num(m,:) = r - lf;

    r = complex(K.Rloss.Q(i),K.Xloss.Q(i));
    lf = (S_Loss-S_Loss2)/imag(S_star(i(1))-S_star2(i(1)));
    data_real_Loss.num(m,:) = [sqrt(mean((real(r) - real(lf)).^2)), max(abs(max((real(r) - real(lf)))), abs(min((real(r) - real(lf))))) ];
    data_imag_Loss.num(m,:) = [sqrt(mean((imag(r) - imag(lf)).^2)), max(abs(max((imag(r) - imag(lf)))), abs(min((imag(r) - imag(lf))))) ];
    data_all_Loss.num(m,:) = r - lf;

case 'Q9'
    i = 9;
    
    r = K.Eabs.Q(:,i);
    lf = (abs(E_star)-abs(E_star2))/imag(mean(S_star(i(:))-S_star2(i(:))));
    data_real_Vabs.num(m,:) = [sqrt(mean((real(r) - real(lf)).^2)), max(abs(max((real(r) - real(lf)))), abs(min((real(r) - real(lf))))) ];
    data_imag_Vabs.num(m,:) = [sqrt(mean((imag(r) - imag(lf)).^2)), max(abs(max((imag(r) - imag(lf)))), abs(min((imag(r) - imag(lf))))) ];
    data_all_Vabs.num(m,:) = r - lf;
    
    r = K.Eang.Q(:,i);
    lf = (angle(E_star)-angle(E_star2))/imag(mean(S_star(i(:))-S_star2(i(:))));
    data_real_Vang.num(m,:) = [sqrt(mean((real(r) - real(lf)).^2)), max(abs(max((real(r) - real(lf)))), abs(min((real(r) - real(lf))))) ];
    data_imag_Vang.num(m,:) = [sqrt(mean((imag(r) - imag(lf)).^2)), max(abs(max((imag(r) - imag(lf)))), abs(min((imag(r) - imag(lf))))) ];
    data_all_Vang.num(m,:) = r - lf;
    
    r = K.Iabs.Q(:,i);
    lf = (abs(I_star)-abs(I_star2))/imag(mean(S_star(i(:))-S_star2(i(:))));
    data_real_I.num(m,:) = [sqrt(mean((real(r) - real(lf)).^2)), max(abs(max((real(r) - real(lf)))), abs(min((real(r) - real(lf))))) ];
    data_imag_I.num(m,:) = [sqrt(mean((imag(r) - imag(lf)).^2)), max(abs(max((imag(r) - imag(lf)))), abs(min((imag(r) - imag(lf))))) ];
    data_all_I.num(m,:) = r - lf;

    r = complex(K.Rloss.Q(i),K.Xloss.Q(i));
    lf = (S_Loss-S_Loss2)/imag(mean(S_star(i(:))-S_star2(i(:))));
    data_real_Loss.num(m,:) = [sqrt(mean((real(r) - real(lf)).^2)), max(abs(max((real(r) - real(lf)))), abs(min((real(r) - real(lf))))) ];
    data_imag_Loss.num(m,:) = [sqrt(mean((imag(r) - imag(lf)).^2)), max(abs(max((imag(r) - imag(lf)))), abs(min((imag(r) - imag(lf))))) ];
    data_all_Loss.num(m,:) = r - lf;
    
case 'P9'
    i = 9;
    r = K.Eabs.P(:,i);
    lf = (abs(E_star)-abs(E_star2))/real(mean(S_star(i(:))-S_star2(i(:))));
    data_real_Vabs.num(m,:) = [sqrt(mean((real(r) - real(lf)).^2)), max(abs(max((real(r) - real(lf)))), abs(min((real(r) - real(lf))))) ];
    data_imag_Vabs.num(m,:) = [sqrt(mean((imag(r) - imag(lf)).^2)), max(abs(max((imag(r) - imag(lf)))), abs(min((imag(r) - imag(lf))))) ];
    data_all_Vabs.num(m,:) = r - lf;
  
    r = K.Eang.P(:,i);
    lf = (angle(E_star)-angle(E_star2))/real(mean(S_star(i(:))-S_star2(i(:))));
    data_real_Vang.num(m,:) = [sqrt(mean((real(r) - real(lf)).^2)), max(abs(max((real(r) - real(lf)))), abs(min((real(r) - real(lf))))) ];
    data_imag_Vang.num(m,:) = [sqrt(mean((imag(r) - imag(lf)).^2)), max(abs(max((imag(r) - imag(lf)))), abs(min((imag(r) - imag(lf))))) ];
    data_all_Vang.num(m,:) = r - lf;
    
    r = K.Iabs.P(:,i);
    lf = (abs(I_star)-abs(I_star2))/real(mean(S_star(i(:))-S_star2(i(:))));
    data_real_I.num(m,:) = [sqrt(mean((real(r) - real(lf)).^2)), max(abs(max((real(r) - real(lf)))), abs(min((real(r) - real(lf))))) ];
    data_imag_I.num(m,:) = [sqrt(mean((imag(r) - imag(lf)).^2)), max(abs(max((imag(r) - imag(lf)))), abs(min((imag(r) - imag(lf))))) ];
    data_all_I.num(m,:) = r - lf;

    r = complex(K.Rloss.P(i),K.Xloss.P(i));
    lf = (S_Loss-S_Loss2)/real(mean(S_star(i(:))-S_star2(i(:))));
    data_real_Loss.num(m,:) = [sqrt(mean((real(r) - real(lf)).^2)), max(abs(max((real(r) - real(lf)))), abs(min((real(r) - real(lf))))) ];
    data_imag_Loss.num(m,:) = [sqrt(mean((imag(r) - imag(lf)).^2)), max(abs(max((imag(r) - imag(lf)))), abs(min((imag(r) - imag(lf))))) ];
    data_all_Loss.num(m,:) = r - lf;
        
case 'E19_abs'
    i = 19;
    
    r = K.Eabs.Eabs(:,i);
    lf = (abs(E_star)-abs(E_star2))/(abs(E_star(i))-abs(E_star2(i)));
    data_real_Vabs.num(m,:) = [sqrt(mean((real(r) - real(lf)).^2)), max(abs(max((real(r) - real(lf)))), abs(min((real(r) - real(lf))))) ];
    data_imag_Vabs.num(m,:) = [sqrt(mean((imag(r) - imag(lf)).^2)), max(abs(max((imag(r) - imag(lf)))), abs(min((imag(r) - imag(lf))))) ];
    data_all_Vabs.num(m,:) = r - lf;
%     
%     r = K.E.Eabs(:,i);
%     lf = (E_star-E_star2)/(abs(E_star(i))-abs(E_star2(i)));
    
    r = K.Eang.Eabs(:,i);
    lf = (angle(E_star)-angle(E_star2))/(abs(E_star(i))-abs(E_star2(i)));
    data_real_Vang.num(m,:) = [sqrt(mean((real(r) - real(lf)).^2)), max(abs(max((real(r) - real(lf)))), abs(min((real(r) - real(lf))))) ];
    data_imag_Vang.num(m,:) = [sqrt(mean((imag(r) - imag(lf)).^2)), max(abs(max((imag(r) - imag(lf)))), abs(min((imag(r) - imag(lf))))) ];
    data_all_Vang.num(m,:) = r - lf;
    
    r = K.Iabs.Eabs(:,i);
    lf = (abs(I_star)-abs(I_star2))/(abs(E_star(i))-abs(E_star2(i)));
    data_real_I.num(m,:) = [sqrt(mean((real(r) - real(lf)).^2)), max(abs(max((real(r) - real(lf)))), abs(min((real(r) - real(lf))))) ];
    data_imag_I.num(m,:) = [sqrt(mean((imag(r) - imag(lf)).^2)), max(abs(max((imag(r) - imag(lf)))), abs(min((imag(r) - imag(lf))))) ];
    data_all_I.num(m,:) = r - lf;

    r = complex(K.Rloss.Eabs(i),K.Xloss.Eabs(i));
    lf = (S_Loss-S_Loss2)/(abs(E_star(i))-abs(E_star2(i)));
    data_real_Loss.num(m,:) = [sqrt(mean((real(r) - real(lf)).^2)), max(abs(max((real(r) - real(lf)))), abs(min((real(r) - real(lf))))) ];
    data_imag_Loss.num(m,:) = [sqrt(mean((imag(r) - imag(lf)).^2)), max(abs(max((imag(r) - imag(lf)))), abs(min((imag(r) - imag(lf))))) ];
    data_all_Loss.num(m,:) = r - lf;

case 'E19_angle'
    i = 19;
    
    r = K.Eabs.Eang(:,i);
    lf = (abs(E_star)-abs(E_star2))/(angle(E_star(i))-angle(E_star2(i)));
    data_real_Vabs.num(m,:) = [sqrt(mean((real(r) - real(lf)).^2)), max(abs(max((real(r) - real(lf)))), abs(min((real(r) - real(lf))))) ];
    data_imag_Vabs.num(m,:) = [sqrt(mean((imag(r) - imag(lf)).^2)), max(abs(max((imag(r) - imag(lf)))), abs(min((imag(r) - imag(lf))))) ];
    data_all_Vabs.num(m,:) = r - lf;
    
    r = K.Eang.Eang(:,i);
    lf = (angle(E_star)-angle(E_star2))/(angle(E_star(i))-angle(E_star2(i)));
    data_real_Vang.num(m,:) = [sqrt(mean((real(r) - real(lf)).^2)), max(abs(max((real(r) - real(lf)))), abs(min((real(r) - real(lf))))) ];
    data_imag_Vang.num(m,:) = [sqrt(mean((imag(r) - imag(lf)).^2)), max(abs(max((imag(r) - imag(lf)))), abs(min((imag(r) - imag(lf))))) ];
    data_all_Vang.num(m,:) = r - lf;
    
    r = K.Iabs.Eang(:,i);
    lf = (abs(I_star)-abs(I_star2))/(angle(E_star(i))-angle(E_star2(i)));
    data_real_I.num(m,:) = [sqrt(mean((real(r) - real(lf)).^2)), max(abs(max((real(r) - real(lf)))), abs(min((real(r) - real(lf))))) ];
    data_imag_I.num(m,:) = [sqrt(mean((imag(r) - imag(lf)).^2)), max(abs(max((imag(r) - imag(lf)))), abs(min((imag(r) - imag(lf))))) ];
    data_all_I.num(m,:) = r - lf;

    r = complex(K.Rloss.Eang(i),K.Xloss.Eang(i));
    lf = (S_Loss-S_Loss2)/(angle(E_star(i))-angle(E_star2(i)));
    data_real_Loss.num(m,:) = [sqrt(mean((real(r) - real(lf)).^2)), max(abs(max((real(r) - real(lf)))), abs(min((real(r) - real(lf))))) ];
    data_imag_Loss.num(m,:) = [sqrt(mean((imag(r) - imag(lf)).^2)), max(abs(max((imag(r) - imag(lf)))), abs(min((imag(r) - imag(lf))))) ];
    data_all_Loss.num(m,:) = r - lf;
        
    
case 'E27'
    i = 27;
    
    r = K.Eabs.Eabs(:,i);
    lf = (abs(E_star)-abs(E_star2))/(E_star(i)-E_star2(i));
    data_real_Vabs.num(m,:) = [sqrt(mean((real(r) - real(lf)).^2)), max(abs(max((real(r) - real(lf)))), abs(min((real(r) - real(lf))))) ];
    data_imag_Vabs.num(m,:) = [sqrt(mean((imag(r) - imag(lf)).^2)), max(abs(max((imag(r) - imag(lf)))), abs(min((imag(r) - imag(lf))))) ];
    data_all_Vabs.num(m,:) = r - lf;
    
    r = K.Eang.Eabs(:,i);
    lf = (angle(E_star)-angle(E_star2))/(E_star(i)-E_star2(i));
    data_real_Vang.num(m,:) = [sqrt(mean((real(r) - real(lf)).^2)), max(abs(max((real(r) - real(lf)))), abs(min((real(r) - real(lf))))) ];
    data_imag_Vang.num(m,:) = [sqrt(mean((imag(r) - imag(lf)).^2)), max(abs(max((imag(r) - imag(lf)))), abs(min((imag(r) - imag(lf))))) ];
    data_all_Vang.num(m,:) = r - lf;
    
    r = K.Iabs.Eabs(:,i);
    lf = (abs(I_star)-abs(I_star2))/(E_star(i)-E_star2(i));
    data_real_I.num(m,:) = [sqrt(mean((real(r) - real(lf)).^2)), max(abs(max((real(r) - real(lf)))), abs(min((real(r) - real(lf))))) ];
    data_imag_I.num(m,:) = [sqrt(mean((imag(r) - imag(lf)).^2)), max(abs(max((imag(r) - imag(lf)))), abs(min((imag(r) - imag(lf))))) ];
    data_all_I.num(m,:) = r - lf;

    r = complex(K.Rloss.Eabs(i),K.Xloss.Eabs(i));
    lf = (S_Loss-S_Loss2)/(E_star(i)-E_star2(i));
    data_real_Loss.num(m,:) = [sqrt(mean((real(r) - real(lf)).^2)), max(abs(max((real(r) - real(lf)))), abs(min((real(r) - real(lf))))) ];
    data_imag_Loss.num(m,:) = [sqrt(mean((imag(r) - imag(lf)).^2)), max(abs(max((imag(r) - imag(lf)))), abs(min((imag(r) - imag(lf))))) ];
    data_all_Loss.num(m,:) = r - lf;
        
end


end

complex(data_real_Vabs.num, data_imag_Vabs.num)
complex(data_real_Vang.num, data_imag_Vang.num)
complex(data_real_I.num, data_imag_I.num)
complex(data_real_Loss.num, data_imag_Loss.num)


%% step 2


[E_star, S_star, I_star] = get_fake_states(idx, Grid_para,Filter_para);
Pslack_1 = real(S_star(1));
Qslack_1 = imag(S_star(1));
PLoss_1 = real(sum(S_star));
QLoss_1 = imag(sum(S_star));
P_1 = real(S_star);
Q_1 = imag(S_star);
Eabs_1 = abs(E_star);
Eang_1 = angle(E_star);
Iabs_1 = abs(I_star);

[E_star, S_star, I_star] = get_fake_states2(idx, Grid_para,Filter_para);
Pslack_2 = real(S_star(1));
Qslack_2 = imag(S_star(1));
PLoss_2 = real(sum(S_star));
QLoss_2 = imag(sum(S_star));
P_2 = real(S_star);
Q_2 = imag(S_star);
Eabs_2 = abs(E_star);
Eang_2 = angle(E_star);
Iabs_2 = abs(I_star);


Eabs_D = Eabs_2 - Eabs_1;
Eabs_Drec =        - K.Eabs.P*(P_1 - P_2)  ...
                   - K.Eabs.Q*(Q_1 - Q_2)  ...
                   - K.Eabs.Eabs*(Eabs_1  - Eabs_2)  ...
                   - K.Eabs.Eang*(Eang_1  - Eang_2); 
max(abs(Eabs_D - Eabs_Drec))


Eang_D = Eang_2 - Eang_1;
Eang_Drec =        - K.Eang.P*(P_1 - P_2)  ...
                   - K.Eang.Q*(Q_1 - Q_2)  ...
                   - K.Eang.Eabs*(Eabs_1  - Eabs_2)  ...
                   - K.Eang.Eang*(Eang_1  - Eang_2); 
max(abs(Eang_D - Eang_Drec))


Iabs_D = Iabs_2 - Iabs_1;
Iabs_Drec =        - K.Iabs.P*(P_1 - P_2)  ...
                   - K.Iabs.Q*(Q_1 - Q_2)  ...
                   - K.Iabs.Eabs*(Eabs_1  - Eabs_2)  ...
                   - K.Iabs.Eang*(Eang_1  - Eang_2); 
max(abs(Iabs_D - Iabs_Drec))

PLoss_D = PLoss_2 - PLoss_1;
PLoss_Drec =       - K.Rloss.P*(P_1 - P_2)  ...
                   - K.Rloss.Q*(Q_1 - Q_2)  ...
                   - K.Rloss.Eabs*(Eabs_1  - Eabs_2)  ...
                   - K.Rloss.Eang*(Eang_1  - Eang_2); 
max(abs(PLoss_D - PLoss_Drec))

QLoss_D = QLoss_2 - QLoss_1;
QLoss_Drec =       - K.Xloss.P*(P_1 - P_2)  ...
                   - K.Xloss.Q*(Q_1 - Q_2)  ...
                   - K.Xloss.Eabs*(Eabs_1  - Eabs_2)  ...
                   - K.Xloss.Eang*(Eang_1  - Eang_2); 
max(abs(QLoss_D - QLoss_Drec))


Pslack_2rec =  - sum(P_2(2:end)) + (PLoss_1 + PLoss_Drec);
mean(abs(Pslack_2 - Pslack_2rec))

Qslack_2rec =  - sum(Q_2(2:end)) + (QLoss_1 + QLoss_Drec);
mean(abs(Qslack_2 - Qslack_2rec))
