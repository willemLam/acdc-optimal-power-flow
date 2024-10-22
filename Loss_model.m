close all
clear all
addpath(genpath(pwd))

%% calibration
afe2_Pref_Pdc = [
-10000  -8200;
-9000   -7440;
-8000   -6680;
-7000   -5900;
-6000   -5000;
-5000   -4240;
-4000   -3540;
-3000   -2730;
-2000   -1910;
-1000   -980;
0	0;
1000	740;
2000	1440;
3000	2000;
4000	2770;
5000	3720;
6000	4620;
7000	5490;
8000	6400;
9000	7250;
10000	8120];

a_dc = afe2_Pref_Pdc(:,1);
b_dc1 = afe2_Pref_Pdc(:,2);
figure
hold on
scatter(afe2_Pref_Pdc(:,1),afe2_Pref_Pdc(:,2))
plot(-10000:10000,0.8119*(-10000:10000) - 193.8)


%%



% initialize plots
set(groot,'defaultFigureVisible','on')
set(0, 'DefaultTextInterpreter', 'Latex')
set(0, 'DefaultLegendInterpreter', 'Latex')
set(0, 'DefaultAxesTickLabelInterpreter', 'Latex')


m_cmap1 = colorgrad(5,'blue_down');
m_cmap2 = colorgrad(5,'red_down');
m_cmap3 = colorgrad(5,'orange_down');

% Initialize grid parameters
tol = 1e-7; n_max = 100;
[Grid_para, Filter_para, idx, idx3, constraints] = initialize_RT("CONNECTED");
Grid_para.Vall_b = [ Grid_para.V_b*ones(Grid_para.n_ac,1) ; Grid_para.Vdc_b*ones(Grid_para.n_dc,1)];
Grid_para.Iall_b = [ Grid_para.I_b*ones(21,1) ; Grid_para.Idc_b*ones(6,1)];

%get the data
FilePath = 'C:\Users\admin\Documents\Willem\acdc_centralized_controller\results\';
FileName_island = fullfile([FilePath, '08-30-2024 16-03_EXPERIMENT_island']);
FileName_connected = fullfile([FilePath, '08-30-2024 16-03_EXPERIMENT_connected']);


solution = Get_all_data_old(FileName_connected);

%% IC1
range = [ 500: 3140 ];
Vdc_ic1 = solution.variables_E(:,22); %fillmissing(Voltage.B23V_mag,'linear')/Grid_para.Vdc_b;
Idc_ic1 = solution.variables_I(:,21); %fillmissing(Current.B23IA_mag,'linear')/Grid_para.Idc_b;
Vac_ic1 = solution.variables_E(:,19); %fillmissing(Voltage.B19V_abs,'linear')/Grid_para.V_b;
Iac_ic1 = solution.variables_I(:,18); %fillmissing(Current.B19IA_mag,'linear')/Grid_para.I_b;
Ploss_ic1 = -(real(solution.variables_S(:,22) + solution.variables_S(:,19)));
Ploss_ic1 = movmean(Ploss_ic1,10);
[a_ac1, b_ac1, c_ac1 , b_dc1, c_dc1] = Losses_parameter_estimation_quadratic_switching(Iac_ic1(range),Idc_ic1(range),Vdc_ic1(range),Ploss_ic1(range),Grid_para);

Ploss_ic1_rec = a_ac1 *abs(Iac_ic1)  +  b_ac1 *abs(Iac_ic1).^2 + c_ac1*Vdc_ic1 + b_dc1*abs(Idc_ic1) + c_dc1 *abs(Idc_ic1).^2;
[a_ac1, b_ac1, c_ac1 , b_dc1, c_dc1]

rmse(Ploss_ic1(range),Ploss_ic1_rec(range))

figure
hold on
plot(abs(Ploss_ic1(range)))
plot(abs(Ploss_ic1_rec(range) ))
legend('observed losses IC1','reconstructed losses')

%% IC2
range = [ 50: 3140 ];
Vdc_ic2 = solution.variables_E(:,23); %fillmissing(Voltage.B23V_mag,'linear')/Grid_para.Vdc_b;
Idc_ic2 = solution.variables_I(:,22); %fillmissing(Current.B23IA_mag,'linear')/Grid_para.Idc_b;
Vac_ic2 = solution.variables_E(:,20); %fillmissing(Voltage.B19V_abs,'linear')/Grid_para.V_b;
Iac_ic2 = solution.variables_I(:,19); %fillmissing(Current.B19IA_mag,'linear')/Grid_para.I_b;
Ploss_ic2 = (real(solution.variables_S(:,22) + solution.variables_S(:,19)));
Ploss_ic2 = movmean(Ploss_ic2,10);
[a_ac2, b_ac2, c_ac2, b_dc2, c_dc2] = Losses_parameter_estimation_quadratic_switching(Idc_ic2(range),Iac_ic2(range),Vdc_ic2(range),Ploss_ic2(range),Grid_para);

Ploss_ic2_rec = a_ac2 *abs(Iac_ic2)  +  b_ac2 *abs(Iac_ic2).^2 + c_ac2*Vdc_ic2 + b_dc2*abs(Idc_ic2) + c_dc2 *abs(Idc_ic2).^2;
[a_ac2, b_ac2, c_ac2, b_dc2, c_dc2]

rmse(Ploss_ic2(range),Ploss_ic2_rec(range))

figure
hold on
plot(Ploss_ic2(range))
plot(Ploss_ic2_rec(range) )
legend('observed losses IC2','reconstructed losses')


%% Sevilla paper
v_ceo = 1;
r_c = 4.5e-3;
v_fo = 2;
r_d = 4.5e-3;

PF = abs(real(solution.variables_S(:,19))./abs(solution.variables_S(:,19)));
m = (sqrt(2)*abs(solution.variables_E(:,19)))./real(solution.variables_E(:,22));
y = m.*PF;
I_0 = sqrt(2) * abs(solution.variables_I(:,19))*Grid_para.I_b;

P_igbt  = I_0 * v_ceo .* (1/(2*pi) + y/8) + r_c*I_0.^2 .* (1/8 + y/(3*pi));
P_diode = I_0 * v_fo  .* (1/(2*pi) - y/8) + r_d*I_0.^2 .* (1/8 - y/(3*pi));

figure
hold on
plot(P_igbt+P_diode)
plot(P_igbt)
plot(P_diode)

figure
plot(m)
%% SE paper
I_0 = abs(solution.variables_I(:,19))*Grid_para.I_b;

P_con_igbt = 3*(v_ceo * 4/pi * I_0 + r_c * I_0.^2);
P_con_diod = 3*(v_fo * 4/pi * I_0 + r_d * I_0.^2);
figure
hold on
plot(P_con_igbt+P_con_diod)
plot(P_con_diod)
plot(P_con_igbt)


Ploss_ic1 = abs(-(real(solution.variables_S(:,22)) + real(solution.variables_S(:,19))))*Grid_para.A_b;
P_sw = Ploss_ic1 - (P_con_igbt+P_con_diod);

range = [ 1: 3100 ];
Idc_ic1 = solution.variables_I(:,18)*Grid_para.I_b;
[a_dc, b_dc, c_dc] = Losses_parameter_estimation(Idc_ic1(range),P_sw(range))

figure
hold on
plot(P_sw)
plot( a_dc* ones(length(Idc_ic1),1) + b_dc* abs(Idc_ic1) + c_dc* abs(Idc_ic1).^2,'LineWidth',2)

figure
scatter(abs(Idc_ic1),P_sw)


figure
hold on
plot(Ploss_ic1)
plot(P_sw)
plot(P_con_igbt+P_con_diod)




%% Visualize data
%power in B19 and B22
figure
Ploss_ic1 = -(real(solution.variables_S(:,22)) + real(solution.variables_S(:,19)));
subplot(2,1,1)
hold on
plot(real(solution.variables_S(:,19)))
plot(real(solution.variables_S(:,22)))
plot(Ploss_ic1)
legend('B19','B22','Delta')
title('power and losses in IC 1')
hold off

subplot(2,1,2)
hold on
plot(abs(solution.variables_I(:,19)))
plot(abs(solution.variables_I(:,22)))
legend('B19','B22')
title('Current in IC 1')
hold off

figure
hold on
scatter(Ploss_ic1,abs(solution.variables_I(:,19)) + abs(solution.variables_I(:,22)))
scatter(Ploss_ic1,abs(solution.variables_I(:,22)))
scatter(Ploss_ic1,abs(solution.variables_I(:,19)))

range = [ 1: 2880 ];
Idc_ic1 = solution.variables_I(:,18);
Iac_ic1 = solution.variables_I(:,21);
Ploss_ic1 = -(abs(real(solution.variables_S(:,22))) - abs(real(solution.variables_S(:,19))));

[a_dc, b_dc1, c_dc1] = Losses_parameter_estimation(Idc_ic1(range),Ploss_ic1(range));
[a_ac1, b_ac1, c_ac1] = Losses_parameter_estimation(Iac_ic1(range),Ploss_ic1(range));
[a_m, b_m, c_m] = Losses_parameter_estimation3(Iac_ic1(range),Idc_ic1(range),Ploss_ic1(range));

figure
hold on
plot(Ploss_ic1*Grid_para.A_b)
% plot( a_dc* ones(length(Idc_ic1),1) + b_dc* abs(Idc_ic1) + c_dc* abs(Idc_ic1).^2,'LineWidth',2)
% plot( a_ac* ones(length(Iac_ic1),1) + b_ac* abs(Iac_ic1) + c_ac* abs(Iac_ic1).^2,'LineWidth',2)
% plot( a* ones(length(Iac_ic1),1) + b* abs(Iac_ic1) + c* abs(Iac_ic1).^2 + d* abs(Idc_ic1) + e* abs(Idc_ic1).^2,'LineWidth',2)
plot( a_m* ones(length(Iac_ic1),1) + b_m* abs(Idc_ic1) + c_m* abs(Iac_ic1).^2,'LineWidth',2)


%%
t0 = '2024-08-30T14:13:00Z'; %-2 for UTC
tf = '2024-08-30T15:39:00Z';
resolution = '1000ms';

Power_idx =   [19,20,23,24];
Voltage_idx = [19,20,23,24];
Current_idx = [19,20,23,24];
Power   = SCADA_Power_Query(Power_idx,t0, tf, resolution,Grid_para);
Voltage   = SCADA_Voltage_Query(Voltage_idx,t0, tf, resolution,Grid_para);
Current = SCADA_Current_Query(Current_idx,t0, tf, resolution,Grid_para);

figure
Ploss_ic1 = - (real(Power.B19P) + real(Power.B23P));
hold on
plot(real(Power.B19P))
plot(real(Power.B23P))
plot(Ploss_ic1)
legend('B19','B22','Delta')
title('power and losses in IC 1')
hold off
%%
%IC1
range = [ 1: 5160 ]; %[ 1250: 3150 ];
Vdc_ic1 = fillmissing(Voltage.B23V_mag,'linear')/Grid_para.Vdc_b;
Idc_ic1 = fillmissing(Current.B23IA_mag,'linear')/Grid_para.Idc_b;
Vac_ic1 = fillmissing(Voltage.B19V_abs,'linear')/Grid_para.V_b;
Iac_ic1 = fillmissing(Current.B19IA_mag,'linear')/Grid_para.I_b;
Ploss_ic1 = ((real(fillmissing(Power.B23P,'linear') + fillmissing(Power.B19P,'linear'))))/Grid_para.A_b;
Ploss_ic1 = movmean(Ploss_ic1,10);
[a_ac1, b_ac1, c_ac1 , b_dc1, c_dc1] = Losses_parameter_estimation_quadratic(Iac_ic1(range),Vac_ic1(range),Idc_ic1(range),Vdc_ic1(range),Ploss_ic1(range));
Ploss_ic1_rec = a_ac1 *abs(Vdc_ic1)  +  b_ac1 *abs(Iac_ic1) + c_ac1 *abs(Iac_ic1).^2 + b_dc1*abs(Idc_ic1) + c_dc1 *abs(Idc_ic1).^2;
[a_ac1, b_ac1, c_ac1 , b_dc1, c_dc1]

rmse(Ploss_ic1,Ploss_ic1_rec)

figure
hold on
plot(Ploss_ic1)
plot(Ploss_ic1_rec )
legend('observed losses IC1','reconstructed losses')


%IC2
% range = [ 1: 5150 ]; %[ 1220: 3250 ];
Vdc_ic2 = fillmissing(Voltage.B24V_mag,'linear')/Grid_para.Vdc_b;
Idc_ic2 = fillmissing(Current.B24IA_mag,'linear')/Grid_para.Idc_b;
Vac_ic2 = fillmissing(Voltage.B20V_abs,'linear')/Grid_para.V_b;
Iac_ic2 = fillmissing(Current.B20IA_mag,'linear')/Grid_para.I_b;
Ploss_ic2 = ((real(fillmissing(Power.B24P,'linear') + fillmissing(Power.B20P,'linear'))))/Grid_para.A_b;
Ploss_ic2 = movmean(Ploss_ic2,1);
[a_ac2, b_ac2, c_ac2 , b_dc2, c_dc2] = Losses_parameter_estimation_quadratic(Iac_ic2(range),Vac_ic2(range),Idc_ic2(range),Vdc_ic2(range),Ploss_ic2(range));
Ploss_ic2_rec = a_ac2 *abs(Vdc_ic2)  +  b_ac2 *abs(Iac_ic2) + c_ac2 *abs(Iac_ic2).^2 + b_dc2*abs(Idc_ic2) + c_dc2 *abs(Idc_ic2).^2;
[a_ac2, b_ac2, c_ac2 , b_dc2, c_dc2]

rmse(Ploss_ic2,Ploss_ic2_rec)

figure
hold on
plot(Ploss_ic2)
plot(Ploss_ic2_rec )
legend('observed losses IC2','reconstructed losses')


[ a_ac1 a_ac2 0 ]
[ b_ac1 b_ac2 0 ]
[ c_ac1 c_ac2 0 ]
[ b_dc1 b_dc2 0 ]
[ c_dc1 c_dc2 0 ]





% [Zlac, Ylac, Zldc, Yldc] = Losses_parameter_estimationpi2([Iac_ic1; Iac_ic2],[Vac_ic1; Vac_ic2],[Idc_ic1; Idc_ic2],[Vdc_ic1; Vdc_ic2],[Ploss_ic1; Ploss_ic2])
% figure
% hold on
% plot([Ploss_ic1; Ploss_ic2])
% plot( Zldc* abs([Idc_ic1; Idc_ic2]).^2 + Yldc* abs([Vdc_ic1; Vdc_ic2]).^2 + Zlac* abs([Iac_ic1; Iac_ic2]).^2 + Ylac* abs([Vac_ic1; Vac_ic2]).^2,'LineWidth',2)

%% opf data
%IC1
range = [10: 1915]; %[ 1220: 3250 ];
Vdc_ic1 = solution.variables_E(:,22);
Idc_ic1 = solution.variables_I(:,21);
Vac_ic1 = solution.variables_E(:,19);
Iac_ic1 = solution.variables_I(:,18);
Ploss_ic1 = -abs((real(solution.variables_S(:,19)) + real(solution.variables_S(:,22))));
Ploss_ic1 = movmean(Ploss_ic1,4);
% [Zlac1, Ylac1, Zldc1, Yldc1] = Losses_parameter_estimationpi2(Iac_ic1(range),Vac_ic1(range),Idc_ic1(range),Vdc_ic1(range),Ploss_ic1(range))
[a, b, c, d, e, f] = Losses_parameter_estimation_c_sw(Iac_ic1(range),Vac_ic1(range),Idc_ic1(range),Vdc_ic1(range),Ploss_ic1(range))
    


figure
hold on
plot(Ploss_ic1)
plot(a*Vdc_ic1 .* abs(Iac_ic1) + b*abs(Iac_ic1) +c* abs(Iac_ic1).^2 + d*Vdc_ic1 .* abs(Idc_ic1) + e*abs(Idc_ic1) + f*abs(Idc_ic1).^2  )
% plot( Zldc1* abs(Idc_ic1).^2 + Yldc1* abs(Vdc_ic1).^2 + Zlac1* abs(Iac_ic1).^2 + Ylac1* abs(Vac_ic1).^2,'LineWidth',2)
legend('observed losses IC1','reconstructed losses')


%IC2
Vdc_ic2 = solution.variables_E(:,23);
Idc_ic2 = solution.variables_I(:,22);
Vac_ic2 = solution.variables_E(:,20);
Iac_ic2 = solution.variables_I(:,19);
Ploss_ic2 = -abs((real(solution.variables_S(:,20)) + real(solution.variables_S(:,23))));
Ploss_ic2 = movmean(Ploss_ic2,1);
[Zlac2, Ylac2, Zldc2, Yldc2] = Losses_parameter_estimationpi2(Iac_ic2(range),Vac_ic2(range),Idc_ic2(range),Vdc_ic2(range),Ploss_ic2(range))

figure
hold on
plot(Ploss_ic2)
plot( Zldc2* abs(Idc_ic2).^2 + Yldc2* abs(Vdc_ic2).^2 + Zlac2* abs(Iac_ic2).^2 + Ylac2* abs(Vac_ic2).^2,'LineWidth',2)
legend('observed losses IC2','reconstructed losses')


figure
hold on
plot(real(solution.variables_I(274:end,[19,20,21]))*Grid_para.A_b)
plot(Power.B19IA_mag(1:2:end))
plot(Power.B20IA_mag(1:2:end))

%% Visualize data
%power in B20 and B23
figure
Ploss_ic2 = -(real(solution.variables_S(:,23)) + real(solution.variables_S(:,20)));
subplot(2,1,1)
hold on
plot(real(solution.variables_S(:,20)))
plot(real(solution.variables_S(:,23)))
plot(Ploss_ic2)
legend('B20','B23','Delta')
title('power and losses in IC 2')
hold off

subplot(2,1,2)
hold on
plot(abs(solution.variables_I(:,20)))
plot(abs(solution.variables_I(:,23)))
legend('B20','B23')
title('Current in IC 2')
hold off

figure
hold on
scatter(Ploss_ic2,abs(solution.variables_I(:,20)) + abs(solution.variables_I(:,23)))
scatter(Ploss_ic2,abs(solution.variables_I(:,23)))
scatter(Ploss_ic2,abs(solution.variables_I(:,20)))
%%
range = [ 815: 1880 ];
Idc_ic2 = solution.variables_I(:,20);
Iac_ic2 = solution.variables_I(:,23);
Ploss_ic1 = -(abs(real(solution.variables_S(:,23))) - abs(real(solution.variables_S(:,20))));

[a_dc, b_dc1, c_dc1] = Losses_parameter_estimation(Idc_ic2(range),Ploss_ic2(range));
[a_ac1, b_ac1, c_ac1] = Losses_parameter_estimation(Iac_ic2(range),Ploss_ic2(range));
[a_m, b_m, c_m] = Losses_parameter_estimation3(Iac_ic2(range),Idc_ic2(range),Ploss_ic2(range));

figure
hold on
plot(Ploss_ic2)
plot( a_m* ones(length(Iac_ic1),1) + b_m* abs(Idc_ic1) + c_m* abs(Iac_ic1).^2,'LineWidth',2)




%% power in B15
for i = 1:length(solution.mode)
    Iflow(:,i) = get_Current_flow(variables.E(:,i),Grid_para);
end
P15 = variables.E(15,:).*conj(Iflow(18,:));
P16 = variables.E(16,:).*conj(Iflow(19,:));
P17 = variables.E(17,:).*conj(Iflow(20,:));

plot(real(P15))
plot(real(P15) + real(variables.S(22,:)))

%% reconstruct setpoint
figure
hold on
plot(real((variables.E(15,:)).*conj(Iflow(19,:)) ))
plot(real((variables.E(19,:)).*conj(Iflow(19,:)) ))
plot(real((variables.E(15,:) - variables.E(19,:)).*conj(Iflow(19,:)) ))

figure
hold on
plot(real(variables.S(20,:)))
plot(real(variables.S(23,:)))
plot(real(solution.P(20,:)))
plot(abs(variables.I(20,:)))
plot(abs(variables.I(23,:)))

plot(real(variables.S(20,:)) + real(variables.S(23,:)));
plot(real(variables.S(21,:)) + real(variables.S(24,:)));


%% Losses model
%switching losses
I_sw = 0.006*abs(variables.I(22:24,:))*Grid_para.Idc_b;
Y4 = I_sw./(750*Grid_para.Ydc_b);
P_sw_dc = real(variables.E(22:24,:)).*real(variables.E(22:24,:)).*Y4;

I_sw = 145*0.006*abs(variables.I(19:21,:))*Grid_para.I_b;
Y4 = I_sw./(750*Grid_para.Ydc_b);
P_sw_ac = real(variables.E(22:24,:)).*real(variables.E(22:24,:)).*Y4;

figure
hold on
plot(Grid_para.A_b*P_sw_dc')
plot(Grid_para.A_b*P_sw_ac')

figure
hold on
plot(abs(variables.I(19:21,:)')*Grid_para.I_b)
plot(abs(variables.I(22:24,:)')*Grid_para.Idc_b)


%Conduction losses
Imag_ac = abs(variables.I(19:21,:));
Imag_dc = abs(variables.I(22:24,:));
DIODE_piecewise = [
0	0
0.01	0.533
25	0.966
50	1.25
75	1.433
100	1.55];
V_diode_forward_dc = interp1(DIODE_piecewise(:,1),DIODE_piecewise(:,2),Imag_dc*Grid_para.Idc_b)./Grid_para.Vdc_b; %pu
Zdiode_c_dc = 4/pi*V_diode_forward_dc./Imag_dc;
V_diode_forward_ac = interp1(DIODE_piecewise(:,1),DIODE_piecewise(:,2),Imag_ac*Grid_para.I_b)./Grid_para.V_b; %pu
Zdiode_c_ac = 4/pi*V_diode_forward_ac./Imag_ac;


IGBT_piecewise = [
0	0
0.01	0.533
25	0.966
50	1.25
75	1.433
100	1.55];
V_igbt_forward_dc = interp1(IGBT_piecewise(:,1),IGBT_piecewise(:,2),Imag_dc*Grid_para.Idc_b)./Grid_para.Vdc_b; %pu
Zigbt_c_dc = 4/pi*V_igbt_forward_dc./Imag_dc;
V_igbt_forward_ac = interp1(IGBT_piecewise(:,1),IGBT_piecewise(:,2),Imag_ac*Grid_para.I_b)./Grid_para.V_b; %pu
Zigbt_c_ac = 4/pi*V_igbt_forward_ac./Imag_ac;
        
P_c_dc = real((variables.I(22:24,:).*(537*(Zdiode_c_dc+Zigbt_c_dc))) .* conj(variables.I(22:24,:))); %22:24
P_c_ac = real((variables.I(19:21,:).*(Zdiode_c_ac+Zigbt_c_ac)) .* conj(variables.I(19:21,:)));


figure
hold on
plot(P_c_dc')
plot(P_c_ac')

%filter
Rc = 2;
C2 = 50e-6;
R1 = 66e-3;
L1 = 0.6e-3;
R2 = 33e-3;
L2 = 0.3e-3;
omega = 2*pi*50;

Za = Grid_para.Y_b*(R1 + 1i*omega*L1);
Zb = Grid_para.Y_b*(R2 + 1i*omega*L2);
Zm = Grid_para.Y_b*(1/ ( 1/(1/ (1i*omega*C2)) + 1/(Rc + 1/ (1i*omega*C2)) ));
Ym = 1/Zm;

Vm = variables.E(19:21,:) + Zb*variables.I(19:21,:);


P_filter_m = Vm .* conj(Vm.*Ym);
P_filter_b = variables.I(19:21,:).* conj(variables.I(19:21,:).*Zb);
P_filter_a = Iflow(16:18,:).* conj(Iflow(16:18,:).*Za);
P_filter = real(P_filter_m + P_filter_b + P_filter_a);


%artificial loss
Z_arti = 0.01;
P_arti = variables.E(22:24,:).*variables.E(22:24,:).*Z_arti;
figure
plot(P_arti'*Grid_para.A_b)

%% plot losses
range = 1:390;
figure
subplot(3,1,1)
hold on
plot(abs(real(variables.S(19,range))+real(variables.S(22,range))*Grid_para.A_b))
plot((1*P_sw_ac(3,range)' + 0*P_c_dc(3,range)')*Grid_para.A_b)
% plot((9*P_sw_ac(1,range)' + 400*P_c_dc(1,range)')*Grid_para.A_b)
ylabel('loss IC1 W')
legend('observed loss','reconstructed loss')
hold off

subplot(3,1,2)
hold on
plot(abs(real(variables.S(20,range))+real(variables.S(23,range))*Grid_para.A_b))
plot((0*P_sw_dc(1,range)' + 1*P_c_dc(1,range)')*Grid_para.A_b)
% plot((1*P_sw_ac(2,range)' + 180*P_c_dc(2,range)')*Grid_para.A_b)
ylabel('loss IC2 W')
legend('observed loss','reconstructed loss')
hold off

subplot(3,1,3)
hold on
plot(abs(real(variables.S(21,range))+real(variables.S(24,range))*Grid_para.A_b))
plot((0*P_sw_ac(2,range)' + 1*P_c_dc(2,range)')*Grid_para.A_b)
ylabel('loss IC3 W')
legend('observed loss','reconstructed loss')
hold off

%% currents
range = 1:390;
figure
subplot(3,1,1)
hold on
plot(abs(abs(variables.I(19,range)))*Grid_para.I_b)
plot(abs(real(variables.I(22,range)))*Grid_para.Idc_b)
hold off

subplot(3,1,2)
hold on
plot(abs(abs(variables.I(20,range)))*Grid_para.I_b)
plot(abs(real(variables.I(23,range)))*Grid_para.Idc_b)
hold off

subplot(3,1,3)
hold on
plot(abs(abs(variables.I(21,range)))*Grid_para.I_b)
plot(abs(real(variables.I(24,range)))*Grid_para.Idc_b)
hold off

% %% Voltages
% range = 1:390;
% figure
% subplot(2,1,1)
% hold on
% plot(abs(abs(variables.E(19,range)))*Grid_para.V_b)
% plot(abs(abs(variables.E(20,range)))*Grid_para.V_b)
% plot(abs(abs(variables.E(21,range)))*Grid_para.V_b)
% 
% hold off
% 
% subplot(2,1,2)
% hold on
% plot(abs(real(variables.E(22,range)))*Grid_para.Vdc_b)
% plot(abs(real(variables.E(23,range)))*Grid_para.Vdc_b)
% plot(abs(real(variables.E(24,range)))*Grid_para.Vdc_b)
% hold off

%% 
% t0 = '2024-08-26T15:37:34Z'; %-2 for UTC
% tf = '2024-08-26T15:50:52Z';

t0 = '2024-08-30T14:25:00Z'; %-2 for UTC
tf = '2024-08-30T14:39:00Z';
resolution = '3000ms';

Power_idx = [1,3,9,11,19,20,21,22,23,24,25,27];
Voltage_idx = [1,15,19,20,21,22,23,24,25,28,29];
Current_idx = [1,19,20,21,22,23,24,25];
Power   = SCADA_Power_Query(Power_idx,t0, tf, resolution,Grid_para);
Voltage = SCADA_Voltage_Query(Voltage_idx,t0, tf, resolution,Grid_para);
Current = SCADA_Current_Query(Current_idx,t0, tf, resolution,Grid_para);
Syncro  = SCADA_Syncro_Query([] ,t0, tf, resolution,Grid_para);
    
%% Losses model
%switching losses
I_sw = 0.006*abs([Current.B23IA_mag,Current.B24IA_mag,Current.B25IA_mag]);
Y4 = I_sw./(750);
P_sw_dc = real([Voltage.B23V_mag,Voltage.B24V_mag,Voltage.B25V_mag]).*real([Voltage.B23V_mag,Voltage.B24V_mag,Voltage.B25V_mag]).*Y4;

I_sw = 0.006*abs([Current.B19IA_mag,Current.B20IA_mag,Current.B21IA_mag]);
Y4 = I_sw./(750);
P_sw_ac = real([Voltage.B23V_mag,Voltage.B24V_mag,Voltage.B25V_mag]).*real([Voltage.B23V_mag,Voltage.B24V_mag,Voltage.B25V_mag]).*Y4;

figure
hold on
plot(P_sw_dc)
plot(P_sw_ac)


%Conduction losses
Imag_dc = abs([Current.B23IA_mag,Current.B24IA_mag,Current.B25IA_mag]);
Imag_ac = abs([Current.B19IA_mag,Current.B20IA_mag,Current.B21IA_mag]);
DIODE_piecewise = [
0	0
0.01	0.533
25	0.966
50	1.25
75	1.433
100	1.55];
V_diode_forward_dc = interp1(DIODE_piecewise(:,1),DIODE_piecewise(:,2),Imag_dc); %V
Zdiode_c_dc = 4/pi*V_diode_forward_dc./Imag_dc; %ohm
V_diode_forward_ac = interp1(DIODE_piecewise(:,1),DIODE_piecewise(:,2),Imag_ac); %V
Zdiode_c_ac = 4/pi*V_diode_forward_ac./Imag_ac; %ohm


IGBT_piecewise = [
0	0
0.01	0.533
25	0.966
50	1.25
75	1.433
100	1.55];
V_igbt_forward_dc = interp1(IGBT_piecewise(:,1),IGBT_piecewise(:,2),Imag_dc); %pu
Zigbt_c_dc = 4/pi*V_igbt_forward_dc./Imag_dc; %ohm
V_igbt_forward_ac = interp1(IGBT_piecewise(:,1),IGBT_piecewise(:,2),Imag_ac); %pu
Zigbt_c_ac = 4/pi*V_igbt_forward_ac./Imag_ac; %ohm
        
P_c_dc = real(([Current.B23IA_mag,Current.B24IA_mag,Current.B25IA_mag].*(Zdiode_c_dc+Zigbt_c_dc)) .* conj([Current.B23IA_mag,Current.B24IA_mag,Current.B25IA_mag])); %22:24
P_c_ac = real(([Current.B19IA_mag,Current.B20IA_mag,Current.B21IA_mag].*(Zdiode_c_ac+Zigbt_c_ac)) .* conj([Current.B19IA_mag,Current.B20IA_mag,Current.B21IA_mag]));


figure
hold on
plot(P_c_dc)
plot(P_c_ac)

%% currents
figure
subplot(3,1,1)
hold on
plot(abs(abs(Current.B19IA_mag)))
plot(abs(abs(Current.B23IA_mag)))
hold off

subplot(3,1,2)
hold on
plot(abs((Current.B20IA_mag)))
plot(abs((Current.B24IA_mag)))
hold off

subplot(3,1,3)
hold on
plot(abs((Current.B21IA_mag)))
plot(abs((Current.B25IA_mag)))
hold off

%% plot losses IC1
figure
subplot(2,1,1)
hold on
plot(abs(real(Power.B19P) + real(Power.B23P)))
plot( 10*P_sw_dc(:,1) + 20*P_c_ac(:,1))
% plot((9*P_sw_ac(1,range)' + 400*P_c_dc(1,range)')*Grid_para.A_b)
ylabel('loss IC1 W')
legend('observed loss','reconstructed loss')
hold off
subplot(2,1,2)
hold on
plot(P_sw_dc(:,1))
plot(P_sw_ac(:,1))
plot(P_c_dc(:,1))
plot(P_c_ac(:,1))
legend('switching dc','switching ac','conduction dc','conduction ac')
ylabel('loss IC1 W')
hold off

%% plot losses IC2
figure
subplot(2,1,1)
hold on
plot(abs((real(Power.B20P) + real(Power.B24P))))
plot( 10*P_sw_dc(:,2) + 20*P_c_ac(:,2))
% plot((9*P_sw_ac(1,range)' + 400*P_c_dc(1,range)')*Grid_para.A_b)
ylabel('loss IC2 W')
legend('observed loss','reconstructed loss')
hold off
subplot(2,1,2)
hold on
plot(P_sw_dc(:,2))
plot(P_sw_ac(:,2))
plot(P_c_dc(:,2))
plot(P_c_ac(:,2))
legend('switching dc','switching ac','conduction dc','conduction ac')
ylabel('loss IC2 W')
hold off

%% plot losses
range = 1:390;
figure
subplot(3,1,1)
hold on
plot(abs(real(Power.B19P) + real(Power.B23P)))
plot((P_sw_dc(:,1)+P_c_ac(:,1)))
% plot((9*P_sw_ac(1,range)' + 400*P_c_dc(1,range)')*Grid_para.A_b)
ylabel('loss IC1 W')
legend('observed loss','reconstructed loss')
hold off

subplot(3,1,2)
hold on
plot(abs(real(Power.B20P)+real(Power.B24P)))
plot((P_sw_dc(:,2)+P_c_ac(:,2)))
% plot((1*P_sw_ac(2,range)' + 180*P_c_dc(2,range)')*Grid_para.A_b)
ylabel('loss IC2 W')
legend('observed loss','reconstructed loss')
hold off

subplot(3,1,3)
hold on
plot(abs(real(Power.B21P)+real(Power.B25P)))
plot((25*P_sw_dc(:,3)+5*P_c_ac(:,3)))
ylabel('loss IC3 W')
legend('observed loss','reconstructed loss')
hold off

%%