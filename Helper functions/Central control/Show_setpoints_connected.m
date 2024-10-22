function Show_setpoints_connected(solution,variables,Grid_para,idx)

%IC
M = [round( solution.Eabs([24:27])' *Grid_para.Vdc_b) ; round( solution.Q([19:22])' *Grid_para.A_b) ; round( solution.P([19:22])' *Grid_para.A_b)];
T1 = array2table(M,'VariableNames',{'IC1','IC2','IC3','IC4'},'RowName',{'DC voltage','Reactive power', 'Active power'}); 


%PV
PVmax =    [solution.PVmax_P_perun;
            solution.PVmax_P_solarmax;
            solution.PVmax_P_facade;
            solution.PVmax_P_emul4;
            solution.PVmax_P_emul5;
            solution.PVmax_P_emul6;
            solution.PVmax_P_emul7];

PV =       [solution.PV_P_perun;
            solution.PV_P_solarmax;
            solution.PVmax_P_facade;
            solution.PV_P_emul4;
            solution.PV_P_emul5;
            solution.PV_P_emul6;
            solution.PV_P_emul7];

M = [PVmax'*Grid_para.A_b; PV'*Grid_para.A_b];
T2 = array2table(M,'VariableNames',{'Perun','Solarmax','Facade','PVemul4','PVemul5','PVemul6','PVemul7'},'RowName',{'MPP','Curtail'}); 


%IC
M = [solution.Iabs(22:28)' *Grid_para.Idc_b ];
T3 = array2table(M,'VariableNames',{'(23 - 27)','(24 - 28)','(25 - 29)','(26 - 30)','(27 - 28)','(28 - 29)','(29 - 30)'},'RowName',{'Current'}); 

%Samsung battery
M = [real(solution.Sambat_P)*Grid_para.A_b real(solution.P(27))*Grid_para.A_b real(variables.P(3))*Grid_para.A_b real(solution.P(1))*Grid_para.A_b; imag(variables.Q(9)) *Grid_para.A_b 0 real(variables.Q(3))*Grid_para.A_b imag(solution.P(1))*Grid_para.A_b];
T4 = array2table(M,'VariableNames',{'Samsung','SupCap','Load','PCC'},'RowName',{'Active power','Reactive power'}); 


disp(T1),disp(T2), disp(T3), disp(T4)

end

