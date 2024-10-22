

SC_soc = 1:100;

for i = 1:100
    Delta_SOC = abs(50 - SC_soc(i));
    SupCap_Pmin_adaptive(i) = max( min(-Delta_SOC * 50/Grid_para.A_b,(-Delta_SOC + 5) * 150/Grid_para.A_b),constraints.SupCap_Pmin);
    SupCap_Pmax_adaptive(i) = min( max( Delta_SOC * 50/Grid_para.A_b,( Delta_SOC - 5) * 150/Grid_para.A_b),constraints.SupCap_Pmax);
end

figure
hold on
plot(SupCap_Pmin_adaptive)
plot(SupCap_Pmax_adaptive)
