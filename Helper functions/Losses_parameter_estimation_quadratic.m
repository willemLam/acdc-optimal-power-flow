function [a, b, c, e ,f] = Losses_parameter_estimation_quadratic(Iac_ic1,Vac_ic1,Idc_ic1,Vdc_ic1,Ploss_ic1,Grid_para)

    H = [ ones(length(Iac_ic1),1), abs(Iac_ic1) , abs(Iac_ic1).^2,  abs(Idc_ic1) , abs(Idc_ic1).^2 ];
    H = [ Vdc_ic1, abs(Iac_ic1).^2, abs(Idc_ic1) , abs(Idc_ic1).^2];
    z = Ploss_ic1;
    R = eye(length(Iac_ic1));

    G = H.' * (R \ H);
    x = G \ (H.' * (R \ z));

    a = x(1);
    b = 0;%x(2);
    c = x(2);
    e = x(3);
    f = x(4);

end