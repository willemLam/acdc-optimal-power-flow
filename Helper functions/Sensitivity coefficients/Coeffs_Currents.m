% REDO EASIER 3 for loops!

function [IKp_abs, IKp, IKq_abs,IKq, IKvabs_abs, IKvabs, IKvang_abs, IKvang, Currents]=Coeffs_Currents(E,SC,Grid_para)


YYL = Grid_para.YYL;
YYT = Grid_para.YYT;
nph = Grid_para.n_ph;
I_indexu = Grid_para.I_indexu;

    Currents = cell(size(YYL)/nph);
    IKp_abs_long = cell(1,size(SC.P,2));
    IKq_abs_long = cell(1,size(SC.P,2));
    IKvabs_abs_long = cell(1,size(SC.P,2));
    IKv_ang_long = cell(1,size(SC.P,2));
    IKp_long = cell(1,size(SC.P,2));
    IKq_long = cell(1,size(SC.P,2));
    IKvabs_long = cell(1,size(SC.P,2));
    IKvang_long = cell(1,size(SC.P,2));
    
    % For all active/reactive power injections
    for l = 1:size(SC.P,2)
        % Re-initialize Temp Matrices
        temp_p = cell(size(YYL)/nph);
        temp_q = cell(size(YYL)/nph); 
        temp_vabs = cell(size(YYL)/nph); 
        temp_vang = cell(size(YYL)/nph); 
        temp_p_complex = cell(size(YYL)/nph);
        temp_q_complex = cell(size(YYL)/nph);
        temp_vabs_complex = cell(size(YYL)/nph);
        temp_vang_complex = cell(size(YYL)/nph);
        
        % For all rows of Admittance Matrix
        for i = 1:nph:size(YYL,1)
            for j = 1:nph:size(YYL,2)
                
                % Compute Currents only once
                if (l == 1)
                    Currents{ceil(i/nph),ceil(j/nph)} = YYL(i:i+(nph-1),j:j+(nph-1))*(E(i:i+(nph-1)) - E(j:j+(nph-1))) + ...
                                    YYT(i:i+(nph-1),j:j+(nph-1))*E(i:i+(nph-1));
                end
                
                if( max(abs(Currents{ceil(i/nph),ceil(j/nph)})) ~= 0)
                    % Compute Complex Coeffs
                    dIij_dPl = YYL(i:i+(nph-1),j:j+(nph-1))*(SC.P(i:i+(nph-1),l) - SC.P(j:j+(nph-1),l)) + ...
                               YYT(i:i+(nph-1),j:j+(nph-1))* SC.P(i:i+(nph-1),l);
                    dIij_dQl = YYL(i:i+(nph-1),j:j+(nph-1))*(SC.Q(i:i+(nph-1),l) - SC.Q(j:j+(nph-1),l)) + ...
                               YYT(i:i+(nph-1),j:j+(nph-1))* SC.Q(i:i+(nph-1),l);
                    dIij_dVabsl = YYL(i:i+(nph-1),j:j+(nph-1))*(SC.Eabs(i:i+(nph-1),l) - SC.Eabs(j:j+(nph-1),l)) + ...
                                  YYT(i:i+(nph-1),j:j+(nph-1))* SC.Eabs(i:i+(nph-1),l);
                    dIij_dVangl = YYL(i:i+(nph-1),j:j+(nph-1))*(SC.Eang(i:i+(nph-1),l) - SC.Eang(j:j+(nph-1),l)) + ...
                                  YYT(i:i+(nph-1),j:j+(nph-1))* SC.Eang(i:i+(nph-1),l);
                           
                    % Compute Coeffs
                    temp_p{ceil(i/nph),ceil(j/nph)} =  ((1./abs(Currents{ceil(i/nph),ceil(j/nph)})).*real(conj(Currents{ceil(i/nph),ceil(j/nph)}).*dIij_dPl));
                    temp_q{ceil(i/nph),ceil(j/nph)} =  ((1./abs(Currents{ceil(i/nph),ceil(j/nph)})).*real(conj(Currents{ceil(i/nph),ceil(j/nph)}).*dIij_dQl));
                    temp_vabs{ceil(i/nph),ceil(j/nph)} =  ((1./abs(Currents{ceil(i/nph),ceil(j/nph)})).*real(conj(Currents{ceil(i/nph),ceil(j/nph)}).*dIij_dVabsl));
                    temp_vang{ceil(i/nph),ceil(j/nph)} =  ((1./abs(Currents{ceil(i/nph),ceil(j/nph)})).*real(conj(Currents{ceil(i/nph),ceil(j/nph)}).*dIij_dVangl));
                    temp_p_complex{ceil(i/nph),ceil(j/nph)} =  dIij_dPl;
                    temp_q_complex{ceil(i/nph),ceil(j/nph)} =  dIij_dQl;
                    temp_vabs_complex{ceil(i/nph),ceil(j/nph)} =  dIij_dVabsl;
                    temp_vang_complex{ceil(i/nph),ceil(j/nph)} =  dIij_dVangl;
                else
                    temp_p{ceil(i/nph),ceil(j/nph)} =  (zeros(nph,1));
                    temp_q{ceil(i/nph),ceil(j/nph)} =  (zeros(nph,1));
                    temp_vabs{ceil(i/nph),ceil(j/nph)} =  (zeros(nph,1));
                    temp_vang{ceil(i/nph),ceil(j/nph)} =  (zeros(nph,1));
                    temp_p_complex{ceil(i/nph),ceil(j/nph)} =  (zeros(nph,1));
                    temp_q_complex{ceil(i/nph),ceil(j/nph)} =  (zeros(nph,1));
                    temp_vabs_complex{ceil(i/nph),ceil(j/nph)} =  (zeros(nph,1));
                    temp_vang_complex{ceil(i/nph),ceil(j/nph)} =  (zeros(nph,1));
                end
                
            end
        end
        
        IKp_abs_long{l} = (temp_p);
        IKq_abs_long{l} = (temp_q);
        IKvabs_abs_long{l} = (temp_vabs)
        IKvang_abs_long{l} = (temp_vang);
        IKp_long{l} = (temp_p_complex);
        IKq_long{l} = (temp_p_complex);
        IKvabs_long{l} = (temp_vabs_complex); 
        IKvang_long{l} = (temp_vang_complex); 

        IKp_abs(:,l) = cell2mat(IKp_abs_long{l}(I_indexu));
        IKq_abs(:,l) = cell2mat(IKq_abs_long{l}(I_indexu));
        IKvabs_abs(:,l) = cell2mat(IKvabs_abs_long{l}(I_indexu));
        IKvang_abs(:,l) = cell2mat(IKvang_abs_long{l}(I_indexu));
        IKp(:,l) = cell2mat(IKp_long{l}(I_indexu));
        IKq(:,l) = cell2mat(IKq_long{l}(I_indexu));
        IKvabs(:,l) = cell2mat(IKvabs_long{l}(I_indexu)); 
        IKvang(:,l) = cell2mat(IKvang_long{l}(I_indexu)); 


    end
        
end
