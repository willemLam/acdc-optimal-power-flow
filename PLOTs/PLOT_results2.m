%%%%%%%%%%%%%%
%% Get data %%
%%%%%%%%%%%%%%
close all
clear all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize grid
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Base values

%% Initialization
tol = 1e-7; n_max = 100;
[Grid_para, Filter_para, idx1, idx3, constraints] = initialize_RT("CONNECTED");
idx = idx1; %single phase equivalent

Grid_para.Vall_b = [ Grid_para.V_b*ones(Grid_para.n_ac,1) ; Grid_para.Vdc_b*ones(Grid_para.n_dc,1)];
Grid_para.Iall_b = [ Grid_para.I_b*ones(21,1) ; Grid_para.Idc_b*ones(6,1)];


FilePath = 'C:\Users\admin\Documents\Willem\acdc_centralized_controller\results\';
% FileName = fullfile([FilePath, '03-04-2024 18-59_SIMULATION.mat']);%4 iter
% FileName = fullfile([FilePath, '03-04-2024 20-31_SIMULATION.mat']);%2 iter
% FileName = fullfile([FilePath, '03-04-2024 20-42_SIMULATION.mat']);%1 iter

% FileName = fullfile([FilePath, '07-26-2024 17-47_EXPERIMENT_connected.mat']);

%% 3AFE + control
FileName = fullfile([FilePath, '08-07-2024 18-09_EXPERIMENT_island']);
% FileName = fullfile([FilePath, '08-07-2024 18-09_EXPERIMENT_connected']);
range = 1496:1700;

s = whos('-file',FileName);
N = length({s.name});
solution = struct();

for i = 5:N
  
if rem(i,100) == 0 %update on loading the data
    i
end
name = ['data_',num2str(i)];
D = load(FileName, genvarname(name));
d = D.(genvarname(name));

    solution.PVmax_P_perun(:,i) = d.PVmax_P_perun;
    solution.PVmax_P_solarmax(:,i) = d.PVmax_P_solarmax;
    solution.PVmax_P_facade(:,i) = d.PVmax_P_facade;

    solution.PV_P_perun(:,i) = d.PV_P_perun;
    solution.PV_P_solarmax(:,i) = d.PV_P_solarmax;
    solution.PV_P_facade(:,i) = d.PV_P_facade;

    solution.P(:,i) = d.P;
    solution.Q(:,i) = d.Q;
    solution.S(:,i) = d.S;
    solution.E(:,i) = d.E;
    solution.Iabs(:,i) = d.Iabs;

    solution.reconstructed.S(:,i) = d.reconstructed.S;
    solution.reconstructed.E(:,i) = d.reconstructed.E;
    solution.reconstructed.Iabs(:,i) = d.reconstructed.Iabs;

    solution.Ploss(:,i) = d.Ploss;
    solution.Qloss(:,i) = d.Qloss;
    solution.Sambat_P(:,i) = d.Sambat_P;
    solution.SupCap_P(:,i) = d.SupCap_P;
    solution.SupCap_SOC(:,i) = d.SupCap_SOC*100;
    solution.IC_Q(:,i) = d.IC_Q;
    solution.IC_E(:,i) = d.IC_E;


    solution.duration(:,i) = d.duration;

    variables.E(:,i) = d.variables.E;
    variables.S(:,i) = d.variables.S;
    variables.I(:,i) = d.variables.I;
end
   

%%%%%%%%%%%%%%%%
%% make plots %%
%%%%%%%%%%%%%%%%

%% CPU time

figure
cdfplot(solution.duration(range))
title('CPU time')
xlabel('Time [sec]')
ylabel('Instances')





%% validate OPF

figure
subplot(5,1,1)
title('Validation SC - mean error btwn model and LF')
plot(rmse(abs(solution.E(:,range)),abs(solution.reconstructed.E(:,range))))
ylabel('voltage magnitude [p.u.]')
subplot(5,1,2)
plot(rmse(angle(solution.E(:,range)),angle(solution.reconstructed.E(:,range))))
ylabel('voltage angle [p.u.]')
subplot(5,1,3)
plot(rmse(abs(solution.Iabs(:,range)),abs(solution.reconstructed.Iabs(:,range))))
ylabel('current magnitude[p.u.]')
subplot(5,1,4)
plot(solution.Ploss(:,range)-sum(real(solution.reconstructed.S(:,range))))
ylabel('loss P [p.u.]')
subplot(5,1,5)
plot(solution.Qloss(:,range)-sum(imag(solution.reconstructed.S(:,range))))
ylabel('loss Q [p.u.]')
xlabel('time')


%% PV curtailment

figure
subplot(1,3,1)
title('Perun')
hold on
ylabel('Power [kW]')
xlabel('time')
plot(solution.PV_P_perun(range),'Linewidth',2)
plot(solution.PVmax_P_perun(range),'Linewidth',2,'LineStyle','--')

subplot(1,3,2)
title('Solarmax')
hold on
ylabel('Power [kW]')
xlabel('time')
plot(solution.PV_P_solarmax(range),'Linewidth',2)
plot(solution.PVmax_P_solarmax(range),'Linewidth',2,'LineStyle','--')

subplot(1,3,3)
title('Solarmax')
hold on
ylabel('Power [kW]')
xlabel('time')
plot(solution.PV_P_facade(range),'Linewidth',2)
plot(solution.PVmax_P_facade(range),'Linewidth',2,'LineStyle','--')

%% Voltage profiles of ICs
figure
plot(solution.reconstructed.E(22:27,range)'*Grid_para.Vdc_b)
title('Voltage profiles')
ylabel('Voltage [p.u.]')
xlabel('time')
ylim([745,755])

%% OPF solution
figure
sgtitle('OPF solution')
idx_P = [1,9,11,19,20,21];
idx_Q = [1,9,11,19,20,21];

subplot(2,1,1)
plot(real(solution.reconstructed.S(idx_P,range))'*Grid_para.A_b,'LineWidth',2)
title('Active power profiles')
legend('GCP','BESS','PV','IC1','IC2','IC3')
ylabel('Active power [p.u.]')
xlabel('time')
xline([68 136],'-',{'Prep. Island','Resync'})

subplot(2,1,2)
plot(imag(solution.reconstructed.S(idx_Q,range))'*Grid_para.A_b,'LineWidth',2)
title('Reactive power profiles')
legend('GCP','BESS','PV','IC1','IC2','IC3')
ylabel('Reactive power [p.u.]')
xlabel('time')


%% state
figure
sgtitle('SCADA state')
idx_P = [1,9,11,19,20,21];
idx_Q = [1,9,11,19,20,21];

subplot(2,1,1)
plot(real(variables.S(idx_P,range))'*Grid_para.A_b,'LineWidth',2)
title('Active power profiles')
legend('GCP','BESS','PV','IC1','IC2','IC3')
ylabel('Active power [p.u.]')
xlabel('time')
xline([68 98 136],'-',{'Prep. Island','Connected but Forming','Resync'})

subplot(2,1,2)
plot(imag(variables.S(idx_Q,range))'*Grid_para.A_b,'LineWidth',2)
title('Reactive power profiles')
legend('GCP','BESS','PV','IC1','IC2','IC3')
ylabel('Reactive power [p.u.]')
xlabel('time')

%% Angle
figure
plot(abs(solution.SupCap_SOC(:,range))')
%%

figure
hold on
plot(real(solution.reconstructed.S(25,range))'*Grid_para.A_b)
plot(real(variables.S(25,range))'*Grid_para.A_b)
xline([68 98 136]+66,'-',{'Prep. Island','Connected but Forming','Resync'})


%% Reactive power losses
figure

subplot(2,1,1)
plot(real(solution.Ploss(range))'*Grid_para.A_b)
title('Active power profiles')
ylabel('Active power loss[p.u.]')
xlabel('time')

subplot(2,1,2)
plot(real(solution.Qloss(range))'*Grid_para.A_b)
title('Reactive power profiles')
ylabel('Reactive power loss[p.u.]')
xlabel('time')

%%
figure
subplot(2,1,1)
plot(real(solution.SupCap_P(range))'*Grid_para.A_b)
title('Active power Super cap')
ylabel('Active power [p.u.]')
xlabel('time')

subplot(2,1,2)
plot(real(solution.SupCap_SOC(range))')
title('Supercap SOC')
ylabel('SOC [%]')
xlabel('time')

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

%% Ampacity limits
figure
title('Ampacity limits')
hold on
ylabel('Current [p.u.]')
xlabel('time')

Congested_lines = [10,];
plot(abs(variables.I(Congested_lines,range)'),'Linewidth',2)
% plot(real(I_line(Congested_lines,:)'),'Linewidth',2)
plot(real(constraints.I_max(Congested_lines)*ones(length(range))),'Linewidth',2,'LineStyle','--')

hold off




