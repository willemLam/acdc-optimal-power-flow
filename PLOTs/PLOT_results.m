%%%%%%%%%%%%%%
%% Get data %%
%%%%%%%%%%%%%%
clear all
warning('off','all')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize grid
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Base values
A_b = 1e5;
V_b= 400;
Y_b = A_b/V_b^2; 

Vdc_b = 800;
Adc_b = A_b;
Ydc_b = Adc_b/Vdc_b^2; 


% FileName = fullfile(['results\09-14-2023 16-32_EXPERIMENT_YES.mat']);
FileName = fullfile(['C:\Users\admin\Documents\Willem\acdc_centralized_controller\results\06-25-2024 23-56_EXPERIMENT.mat']);

% results\07-Sep-2023_EXPERIMENT_6.mat

% results\08-Sep-2023_EXPERIMENT_6.mat
% results\11-Sep-2023_EXPERIMENT_1.mat
% results\11-Sep-2023_SIMULATION_1.mat

s = whos('-file',FileName);
N = length({s.name});

for i = 1:N
name = ['data_',num2str(i)];
D = load(FileName, genvarname(name));
d= D.(genvarname(name));

%     E_opt_IC(:,i) = d.E_optimal_IC;
%     Q_opt_IC(:,i) = imag(d.S_optimal_IC);
%     P_DCT(:,i) = d.S_optimal_DCT;
%     P_opt_PV(:,i) = d.S_optimal_PV;
%     S_max_PV(:,i) = d.S_max_PV;
    CPU_time(i) = d.duration;
%     E(:,i) = d.E_optimal;
%     S(:,i) = d.S_optimal;
%     I(:,i) = d.I_optimal;
%     E_bus(:,i) = d.E_scada;
%     S_bus(:,i) = d.S_scada;
%     I_line(:,i) = d.I_scada;
%     Amp(:,i) = d.Ampacities;
%     Sp_checkk(i) = d.Sp_check;
%     Spn_DCT(i,:) = d.Spn_optimal_DCT;
% 
%     P_DCT_actual(i,:) = d.P_DCT_actual;
%     P_DCT_expected(i,:) = d.P_DCT_expected; 
%     
%     E_DCT_actual(i,:) = d.E_DCT_actual;
%     E_DCT_expected(i,:) = d.E_DCT_expected; 

end
   
%%%%%%%%%%%%%%%%
%% make plots %%
%%%%%%%%%%%%%%%%

%% CPU time

figure
cdfplot(CPU_time(5:end))
title('CPU time')
xlabel('Time [sec]')
ylabel('Instances')

%% validate OPF

figure
subplot(3,1,1)
title('Validation SC - mean error btwn model and LF')
plot(mean(abs(E-E_bus)))
ylabel('voltage [p.u.]')
subplot(3,1,2)
plot(mean(abs(I-I_line)))
ylabel('current [p.u.]')
subplot(3,1,3)
plot(mean(abs(S-S_bus)))
ylabel('power [p.u.]')
xlabel('time')

% ylim([-0.01,0.01])


%% Curtailment

figure
title('Curtailment')
hold on
ylabel('Power [kW]')
xlabel('time')

plot(real(A_b/1000*P_opt_PV'),'Linewidth',2)
plot(real(A_b/1000*S_max_PV(2:end,:)'),'Linewidth',2,'LineStyle','--')



%% Validate Facade
% figure
% title('Curtailment')
% hold on
% ylabel('Power [kW]')
% xlabel('time')
% 
% plot(real(A_b/1000*real(S_bus(11,:))'-S_max_PV(2,:)'),'Linewidth',2)
% plot(real(A_b/1000*S_max_PV(3,:)'),'Linewidth',2,'LineStyle','--')
% legend('Aggregated Facade','facade max')

%% DCT power flow
figure
title('DCT power flow')
hold on
ylabel('Power [kW]')
xlabel('time')

plot(real(A_b/1000*P_DCT([1,2],:)'),'Linewidth',2)


%% Ampacity limits
figure
title('Ampacity limits')
hold on
ylabel('Current [p.u.]')
xlabel('time')

Congested_lines = [10,];
plot(real(I(Congested_lines,:)'),'Linewidth',2)
plot(real(I_line(Congested_lines,:)'),'Linewidth',2)
plot(real(Amp(Congested_lines,:)'),'Linewidth',2,'LineStyle','--')



%% Voltage profiles
figure
plot(abs(E)')
title('Voltage profiles')
ylabel('Voltage [p.u.]')
xlabel('time')


%% Power profiles
figure
subplot(1,2,1)
plot(real(S)')
title('P star')
ylabel('Active power [p.u.]')
xlabel('time')
subplot(1,2,2)
plot(imag(S)')
title('Q star')
ylabel('Reactive power [p.u.]')
xlabel('time')

%% Validity of the solution
figure
plot(Sp_checkk)
title('Valid solution?: 1-YES   0-NO')
ylabel('binary')
xlabel('time')

%% DCT loss
S_bus_rec = E_bus.*conj(Grid_para.YY*E_bus);
figure
hold on
plot(real(S_bus_rec(24,:)))
plot(-real(S_bus_rec(25,:))+0.0008)
title('DC loss')
ylabel('binary')
xlabel('time')

 
figure
hold on
scatter((E_DCT_actual(:,1)-E_DCT_actual(:,2))*Vdc_b,P_DCT_actual*A_b)
scatter((E_DCT_expected(:,1)-E_DCT_expected(:,2))*Vdc_b,P_DCT_expected*A_b)
title('DCT')
ylabel('P')
xlabel('delta Vdc')

figure
hold on
scatter((E_DCT_actual(:,1)  -E_DCT_actual(:,2))*Vdc_b,  (P_DCT_actual(:,1)+  P_DCT_actual(:,2))*A_b)
scatter((E_DCT_expected(:,1)-E_DCT_expected(:,2))*Vdc_b,(P_DCT_expected(:,1)+P_DCT_expected(:,2))*A_b)
title('DCT')
ylabel('P loss')
xlabel('delta Vdc')

