
function [E_star,S_star,I_star] = get_fake_states2(idx, Grid_para,Filter_para)


E = zeros(Grid_para.n_nodes,1);
S = zeros(Grid_para.n_nodes,1);
E_0 = ones(Grid_para.n_nodes,1);
E(1) = 0.995 + 1i*0.1;
E(23:26) = 0.9377;
E(27) = 0.9377;
E(19) = 0.99 + 1i*0.01;
S(3) = -0.05 + 1i*0.01;
S(9) = 0.192 + 1i*0.006;
S(11) = 0.082 + 1i*0.042;
S(14) = -0.2 - 1i*0.1;
S(28:30) = 0.032;
S(19:22) = 1i*0.092;


tol = 1e-6;
n_max = 20;



[E_star,J,n_iter] = NR_rectangularACDC_1ph_general_V2(Grid_para,Filter_para,S,E,E_0,idx,tol,n_max);

S_star = E_star.*conj(Grid_para.YY*E_star);
I_star = get_Current_flow(E_star,Grid_para);



end