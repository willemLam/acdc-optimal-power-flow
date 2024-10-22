function [Zlac, Ylac,Zldc, Yldc] = Losses_parameter_estimationpi2(Iac_ic1,Vac_ic1,Idc_ic1,Vdc_ic1,Ploss_ic1)
    
    H = [ abs(Iac_ic1).^2 , abs(Vac_ic1).^2 , abs(Idc_ic1).^2 , abs(Vdc_ic1).^2 ];
%     H = [ abs(Iac_ic1).^2 , abs(Vac_ic1).^2 , abs(Idc_ic1).^2 ];
      z = Ploss_ic1;
    R = eye(length(Iac_ic1));

    G = H.' * (R \ H);
    x = G \ (H.' * (R \ z));
    Zlac = x(1);
    Ylac = x(2);
    Zldc = x(3);
    Yldc = x(4);

end