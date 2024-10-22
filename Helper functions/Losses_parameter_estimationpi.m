function [Zl, Yl] = Losses_parameter_estimationpi(Iac_ic1,Vac_ic1,Ploss_ic1)
    
    H = [ abs(Iac_ic1).^2 , abs(Vac_ic1).^2 ];
   
    z = Ploss_ic1;
    R = eye(length(Iac_ic1));

    G = H.' * (R \ H);
    x = G \ (H.' * (R \ z));
    Zl = x(1);
    Yl = x(2);

end