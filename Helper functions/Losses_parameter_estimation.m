function [a ,b, c] = Losses_parameter_estimation(Iac_ic1,Ploss_ic1)
    
    H = [ ones(length(Iac_ic1),1) , abs(Iac_ic1) , abs(Iac_ic1).^2 ];
    z = Ploss_ic1;
    R = eye(length(Iac_ic1));

    G = H.' * (R \ H);
    x = G \ (H.' * (R \ z));
    a = x(1);
    b = x(2);
    c = x(3);

end