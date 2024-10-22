function [C_rp , C_rq, C_rvabs, C_rvang, C_xp, C_xq, C_xvabs, C_xvang] =Coeffs_Losses(E_LF, K_E, Grid_para)

YY = Grid_para.YY;

temp_p = zeros(1,size(K_E.P,2));
temp_q = zeros(1,size(K_E.P,2));
temp_vabs = zeros(1,size(K_E.P,2));
temp_vang = zeros(1,size(K_E.P,2));


for l = 1:size(K_E.P,2)
    
    for i = 1:size(E_LF,1)
        
        tmp_YijEj = 0;
        tmp_YijdEjdP = 0;
        tmp_YijdEjdQ = 0;
        tmp_YijdEjdVr = 0;
        tmp_YijdEjdVi = 0;
        
        for j = 1:size(E_LF,1)
            tmp_YijEj = tmp_YijEj + conj(YY(i,j))*conj(E_LF(j));
            tmp_YijdEjdP = tmp_YijdEjdP + conj(YY(i,j))*conj(K_E.P(j,l));
            tmp_YijdEjdQ = tmp_YijdEjdQ + conj(YY(i,j))*conj(K_E.Q(j,l));
            tmp_YijdEjdVr = tmp_YijdEjdVr + conj(YY(i,j))*conj(K_E.Eabs(j,l));
            tmp_YijdEjdVi = tmp_YijdEjdVi + conj(YY(i,j))*conj(K_E.Eang(j,l));
        end 
         
        temp_p(l) = temp_p(l) + K_E.P(i,l) *tmp_YijEj + E_LF(i)*tmp_YijdEjdP;   
        temp_q(l) = temp_q(l) + K_E.Q(i,l)*tmp_YijEj + E_LF(i)*tmp_YijdEjdQ;
        temp_vabs(l) = temp_vabs(l) + K_E.Eabs(i,l)*tmp_YijEj + E_LF(i)*tmp_YijdEjdVr;
        temp_vang(l) = temp_vang(l) + K_E.Eang(i,l)*tmp_YijEj + E_LF(i)*tmp_YijdEjdVi;
    end 
end

C_rp = real(temp_p); 
C_rq = real(temp_q);
C_rvabs = real(temp_vabs);
C_rvang = real(temp_vang);

C_xp = imag(temp_p);
C_xq = imag(temp_q);
C_xvabs = imag(temp_vabs);
C_xvang = imag(temp_vang);
end