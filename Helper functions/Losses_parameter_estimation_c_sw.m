function [a, b, c, d, e ,f] = Losses_parameter_estimation_c_sw(Iac_ic1,Vac_ic1,Idc_ic1,Vdc_ic1,Ploss_ic1)
    
    H = [ Vdc_ic1 .* abs(Iac_ic1) , abs(Iac_ic1) , abs(Iac_ic1).^2, Vdc_ic1 .* abs(Idc_ic1) , abs(Idc_ic1) , abs(Idc_ic1).^2  ];
    
    H = [   ones(length(Iac_ic1),1), abs(Iac_ic1) , abs(Iac_ic1).^2, ones(length(Idc_ic1),1), abs(Idc_ic1) , abs(Idc_ic1).^2, ];
    H = [   abs(Vac_ic1).^2 , abs(Iac_ic1).^2, abs(Vdc_ic1).^2, abs(Idc_ic1).^2  ];
   

    z = Ploss_ic1;
    R = eye(length(Iac_ic1));

    G = H.' * (R \ H);
    x = G \ (H.' * (R \ z));

    a = x(1);
    b = 0;
    c = x(2);
    d = x(3);
    e = 0;
    f = x(4);

end