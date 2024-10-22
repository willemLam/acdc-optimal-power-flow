function Show_setpoints(solution,variables,constraints,Grid_para,idx)

%IC
M = [ round( real(solution.reconstructed.S([19:21]))' *Grid_para.A_b),0 , 0 ; round( imag(solution.reconstructed.S(19:21))' *Grid_para.A_b), 0 , 0; round(Grid_para.A_b*solution.P_quad_ic,0)' 0 ,0; round( solution.reconstructed.E(22:24)' *Grid_para.Vdc_b,1),0 , 0; round(solution.reconstructed.Iabs(18:20)' *Grid_para.I_b,2),0 , round(abs((solution.reconstructed.Iabs(10))*Grid_para.I_b),2) ; 0,0,0,0,round(constraints.I_max(10)*Grid_para.I_b,2) ];
T1 = array2table(M,'VariableNames',{'IC1','IC2','IC3',' ','B10-B11' },'RowName',{'P', 'Q','Ploss','Vdc','Iac','Imax'}); 

%Samsung battery
% M = round([(solution.Sambat_P) (solution.SupCap_P) (solution.P(1)) (solution.PV_P_solarmax) (solution.PV_P_perun); (solution.Q(9)) 0 (solution.Q(1)) (solution.Q(11)) (solution.Q(9)); 0 100*(solution.SupCap_E-0.5e-3)/1.5e-3/Grid_para.A_b  0 0 0]*Grid_para.A_b);
M = round([(solution.Sambat_P)*Grid_para.A_b (solution.SupCap_P)*Grid_para.A_b (solution.P(1))*Grid_para.A_b (variables.PVmax_P_solarmax)*Grid_para.A_b (solution.PV_P_perun)*Grid_para.A_b (variables.PV_P_facade)*Grid_para.A_b; (solution.Q(9))*Grid_para.A_b 0 (solution.Q(1))*Grid_para.A_b (solution.Q(11))*Grid_para.A_b (solution.Q(9))*Grid_para.A_b 0; 0 100*(solution.SupCap_SOC)  0 0 0 0; 0 solution.reconstructed.E(25)' *Grid_para.Vdc_b 0 0 0 0; 0 0 round(solution.reconstructed.Iabs(1)*Grid_para.I_b,2) 0 0 0 ]);
T2 = array2table(M,'VariableNames',{'BESS','SupCap','PCC','PV_sol','PV_perun','PV_facade'},'RowName',{'P','Q','SOC','Vdc','Iac'}); 





disp(T1),disp(T2)

end

