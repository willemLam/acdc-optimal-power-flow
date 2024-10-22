function [a, b, c, e ,f] = Losses_parameter_estimation_quadratic_switching(Iac_ic1,Idc_ic1,Vdc_ic1,Ploss_ic1,Grid_para)
    
    v_ceo = 0*1/Grid_para.V_b;
    r_c = 4.5e-3*Grid_para.Y_b;

    P_cond = v_ceo * abs(Iac_ic1)  +  r_c * abs(Iac_ic1).^2;

    H = [ ones(length(Idc_ic1),1).*Vdc_ic1,  abs(Idc_ic1) , abs(Idc_ic1).^2 ];


    z = Ploss_ic1 - P_cond;
    R = eye(length(Iac_ic1));

    G = H.' * (R \ H);
    x = G \ (H.' * (R \ z));

    a = v_ceo;
    b = r_c;
    c = x(1);
    e = x(2);
    f = x(3);

end