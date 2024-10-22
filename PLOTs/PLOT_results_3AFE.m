%%%%%%%%%%%%%%
%% Get data %%
%%%%%%%%%%%%%%
close all
clear all
addpath(genpath(pwd))



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
% FileName = fullfile([FilePath, '08-07-2024 18-09_EXPERIMENT_connected']);
% range = 1450:1700;

FileName = fullfile([FilePath, '08-26-2024 17-37_EXPERIMENT_island']);
FileName = fullfile([FilePath, '08-26-2024 17-37_EXPERIMENT_connected']);
range = 130:370;


% FileName_island = fullfile([FilePath, '08-27-2024 10-38_EXPERIMENT_island']);
% FileName_connected = fullfile([FilePath, '08-27-2024 10-38_EXPERIMENT_connected']);
% range = 5:90;

    s = whos('-file',FileName);
    N = length({s.name});
    solution = struct();
    
    for i = 1:N
      
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
    %     solution.IC_E(:,i) = d.IC_E;
        solution.feasible(i) = d.proceed;
    
%         solution.mode{:,i} = d.mode;
%         solution.time{:,i} = d.time;
    
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
hold on
cdfplot(solution.duration(range))
legend('island','connected')
title('CPU time')
xlabel('Time [sec]')
ylabel('Instances')


%% SOC

pks_range = 135:330;
figure
hold on
plot(solution.SupCap_SOC(pks_range))
title('SOC')
xlabel('Time [sec]')
ylabel('SOC ')


[pks,locs] = findpeaks(solution.SupCap_SOC(pks_range),1:length(pks_range),'MinPeakProminence',0.5);
% visualizse the selected peaks
plot(locs,pks,'rv')

% calculate the approximate frequency
meanCycle = (locs(end)-locs(1))/(length(pks)-1); %3.2586 in sec

solution.SupCap_SOC_mvavg = movmean(solution.SupCap_SOC(pks_range),meanCycle);
plot(solution.SupCap_SOC_mvavg)
plot(movmean(solution.SupCap_SOC(pks_range),1))
plot(movmean(solution.SupCap_SOC(pks_range),2))
plot(movmean(solution.SupCap_SOC(pks_range),3))
plot(movmean(solution.SupCap_SOC(pks_range),4))

locs_mean = locs(1):meanCycle:locs(end);
plot(locs_mean,pks,'x')

%% 

% [pks2,locs2] = findpeaks(solution.SupCap_SOC_mvavg,1:length(solution.SupCap_SOC_mvavg),'MinPeakProminence',0.05);
% % visualizse the selected peaks
% plot(locs2,pks2,'rv')
% 
% % calculate the approximate frequency
% meanCycle2 = (locs2(end)-locs2(1))/(length(pks2)); %in sec
% 
% locs_mean2 = locs2(1):meanCycle2:locs2(end);
% plot(locs_mean2,pks,'x')




hold off


%% validate OPF

figure
subplot(5,1,1)
title('Validation SC - mean error btwn model and LF')
plot(rmse(abs(solution.E),abs(solution.reconstructed.E)))
ylabel('voltage magnitude [p.u.]')
subplot(5,1,2)
plot(rmse(angle(solution.E),angle(solution.reconstructed.E)))
ylabel('voltage angle [p.u.]')
subplot(5,1,3)
plot(rmse(abs(solution.Iabs),abs(solution.reconstructed.Iabs)))
ylabel('current magnitude[p.u.]')
subplot(5,1,4)
plot(solution.Ploss-sum(real(solution.reconstructed.S)))
ylabel('loss P [p.u.]')
subplot(5,1,5)
plot(solution.Qloss-sum(imag(solution.reconstructed.S)))
ylabel('loss Q [p.u.]')
xlabel('time')


%% PV curtailment

figure
subplot(1,3,1)
title('Perun')
hold on
ylabel('Power [kW]')
xlabel('time')
plot(solution.PV_P_perun,'Linewidth',2)
plot(solution.PVmax_P_perun,'Linewidth',2,'LineStyle','--')

subplot(1,3,2)
title('Solarmax')
hold on
ylabel('Power [kW]')
xlabel('time')
plot(solution.PV_P_solarmax,'Linewidth',2)
plot(solution.PVmax_P_solarmax,'Linewidth',2,'LineStyle','--')

subplot(1,3,3)
title('Solarmax')
hold on
ylabel('Power [kW]')
xlabel('time')
plot(solution.PV_P_facade,'Linewidth',2)
plot(solution.PVmax_P_facade,'Linewidth',2,'LineStyle','--')

%% Voltage profiles of ICs
% figure
% plot(solution.reconstructed.E(23:30,:)'*Grid_para.Vdc_b)
% title('Voltage profiles')
% ylabel('Voltage [p.u.]')
% xlabel('time')
% ylim([745,755])

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

%% OPF solution feasibiltiy
figure
plot(solution.feasible)
%% Difference in state and OPF

% Define the pastel colors for shading each mode
shading_colors = containers.Map({'Grid connected', 'Prepare Island', 'Forming but connected', 'Island', 'Prepare for connection'}, ...
                                {[0.6 0.7 0.9], [1.0 0.8 0.8], [1.0 0.6 0.6], [0.8 0.4 0.4], [0.7 0.8 1.0]});

figure
sgtitle('Reactive powers - difference State and OPF')
subplot(5,1,1)
hold on 
plot(imag(solution.reconstructed.S(9,range))'*Grid_para.A_b,'LineWidth',2)
plot(imag(variables.S(9,range))'*Grid_para.A_b,'LineWidth',2)
for j = 1:length(range)
    mode = solution.mode{range(j)};
    x = [j j j+1 j+1];    y = ylim;            % X and Y limits of the plot for full height shading
    fill(x, [y(1) y(2) y(2) y(1)], shading_colors(mode), 'EdgeColor', 'none', 'FaceAlpha', 0.9);
end
legend('BESS-OPF','BESS-state')


subplot(5,1,2)
hold on 
plot(imag(solution.reconstructed.S(11,range))'*Grid_para.A_b,'LineWidth',2)
plot(imag(variables.S(11,range))'*Grid_para.A_b,'LineWidth',2)
for j = 1:length(range)
    mode = solution.mode{range(j)};
    x = [j j j+1 j+1];    y = ylim;            % X and Y limits of the plot for full height shading
    fill(x, [y(1) y(2) y(2) y(1)], shading_colors(mode), 'EdgeColor', 'none', 'FaceAlpha', 0.3);
end
legend('PV-OPF','PV-state')

subplot(5,1,3)
hold on 
plot(imag(solution.reconstructed.S(19,range))'*Grid_para.A_b,'LineWidth',2)
plot(imag(variables.S(19,range))'*Grid_para.A_b,'LineWidth',2)
for j = 1:length(range)
    mode = solution.mode{range(j)};
    x = [j j j+1 j+1];    y = ylim;            % X and Y limits of the plot for full height shading
    fill(x, [y(1) y(2) y(2) y(1)], shading_colors(mode), 'EdgeColor', 'none', 'FaceAlpha', 0.3);
end
legend('IC1-OPF','IC1-state')

subplot(5,1,4)
hold on 
plot(imag(solution.reconstructed.S(20,range))'*Grid_para.A_b,'LineWidth',2)
plot(imag(variables.S(20,range))'*Grid_para.A_b,'LineWidth',2)
for j = 1:length(range)
    mode = solution.mode{range(j)};
    x = [j j j+1 j+1];    y = ylim;            % X and Y limits of the plot for full height shading
    fill(x, [y(1) y(2) y(2) y(1)], shading_colors(mode), 'EdgeColor', 'none', 'FaceAlpha', 0.3);
end
legend('IC2-OPF','IC2-state')

subplot(5,1,5)
hold on 
plot(imag(solution.reconstructed.S(21,range))'*Grid_para.A_b,'LineWidth',2)
plot(imag(variables.S(21,range))'*Grid_para.A_b,'LineWidth',2)
for j = 1:length(range)
    mode = solution.mode{range(j)};
    x = [j j j+1 j+1];    y = ylim;            % X and Y limits of the plot for full height shading
    fill(x, [y(1) y(2) y(2) y(1)], shading_colors(mode), 'EdgeColor', 'none', 'FaceAlpha', 0.3);
end
legend('IC3-OPF','IC3-state')


%% Difference in state and OPF

figure

sgtitle('Active powers - difference State and OPF')
subplot(5,1,1)
hold on 
plot(real(solution.reconstructed.S(9,range))'*Grid_para.A_b,'LineWidth',2)
plot(real(variables.S(9,range))'*Grid_para.A_b,'LineWidth',2)
for j = 1:length(range)
    mode = solution.mode{range(j)};
    x = [j j j+1 j+1];    y = ylim;            % X and Y limits of the plot for full height shading
    fill(x, [y(1) y(2) y(2) y(1)], shading_colors(mode), 'EdgeColor', 'none', 'FaceAlpha', 0.3);
end
legend('BESS-OPF','BESS-state')

subplot(5,1,2)
hold on 
plot(real(solution.reconstructed.S(11,range))'*Grid_para.A_b,'LineWidth',2)
plot(real(variables.S(11,range))'*Grid_para.A_b,'LineWidth',2)
for j = 1:length(range)
    mode = solution.mode{range(j)};
    x = [j j j+1 j+1];    y = ylim;            % X and Y limits of the plot for full height shading
    fill(x, [y(1) y(2) y(2) y(1)], shading_colors(mode), 'EdgeColor', 'none', 'FaceAlpha', 0.3);
end
legend('PV-OPF','PV-state')

subplot(5,1,3)
hold on 
plot(real(solution.reconstructed.S(19,range))'*Grid_para.A_b,'LineWidth',2)
plot(real(variables.S(19,range))'*Grid_para.A_b,'LineWidth',2)
for j = 1:length(range)
    mode = solution.mode{range(j)};
    x = [j j j+1 j+1];    y = ylim;            % X and Y limits of the plot for full height shading
    fill(x, [y(1) y(2) y(2) y(1)], shading_colors(mode), 'EdgeColor', 'none', 'FaceAlpha', 0.3);
end
legend('IC1-OPF','IC1-state')

subplot(5,1,4)
hold on 
plot(real(solution.reconstructed.S(20,range))'*Grid_para.A_b,'LineWidth',2)
plot(real(variables.S(20,range))'*Grid_para.A_b,'LineWidth',2)
for j = 1:length(range)
    mode = solution.mode{range(j)};
    x = [j j j+1 j+1];    y = ylim;            % X and Y limits of the plot for full height shading
    fill(x, [y(1) y(2) y(2) y(1)], shading_colors(mode), 'EdgeColor', 'none', 'FaceAlpha', 0.3);
end
legend('IC2-OPF','IC2-state')

subplot(5,1,5)
hold on 
plot(real(solution.reconstructed.S(21,range))'*Grid_para.A_b,'LineWidth',2)
plot(real(variables.S(21,range))'*Grid_para.A_b,'LineWidth',2)
for j = 1:length(range)
    mode = solution.mode{range(j)};
    x = [j j j+1 j+1];    y = ylim;            % X and Y limits of the plot for full height shading
    fill(x, [y(1) y(2) y(2) y(1)], shading_colors(mode), 'EdgeColor', 'none', 'FaceAlpha', 0.3);
end
legend('IC3-OPF','IC3-state')



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

%% SOC
figure
plot(abs(solution.SupCap_SOC(:,range))')
%% Reactive power losses
figure

subplot(2,1,1)
plot(real(solution.Ploss)'*Grid_para.A_b)
title('Active power profiles')
ylabel('Active power loss[p.u.]')
xlabel('time')

subplot(2,1,2)
plot(real(solution.Qloss)'*Grid_para.A_b)
title('Reactive power profiles')
ylabel('Reactive power loss[p.u.]')
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
% figure
% title('Ampacity limits')
% hold on
% ylabel('Current [p.u.]')
% xlabel('time')
% 
% Congested_lines = [10,];
% plot(real(I(Congested_lines,:)'),'Linewidth',2)
% plot(real(I_line(Congested_lines,:)'),'Linewidth',2)
% plot(real(Amp(Congested_lines,:)'),'Linewidth',2,'LineStyle','--')


%% influx

    t0 = '2024-08-07T17:01:48Z'; %-2 for UTC
    tf = '2024-08-07T17:05:59Z';
    resolution = '1000ms';

    Power_idx = [1,3,9,11,19,20,21,27];
    Voltage_idx = [1,15,23,24,28,29];
    Current_idx = [1,19];
    Power   = SCADA_Power_Query(Power_idx,t0, tf, resolution,Grid_para);
    Voltage = SCADA_Voltage_Query(Voltage_idx,t0, tf, resolution,Grid_para);
    Current = SCADA_Current_Query(Current_idx,t0, tf, resolution,Grid_para);
    Syncro  = SCADA_Syncro_Query([] ,t0, tf, resolution,Grid_para);
    
    figure
    hold on    
    plot(Power.Time,Power.B19P)
    plot(Power.Time,Power.B20P)
    plot(Power.Time,Power.B21P)
    legend({'B19P','B20P','B21P'})

    figure
    hold on
    plot(Current.time,Current.B01IA_mag)
    plot(Current.time,Current.B01IB_mag)
    plot(Current.time,Current.B01IC_mag)
    plot(Current.time,Current.B19IA_mag)
    plot(Current.time,Current.B19IB_mag)
    plot(Current.time,Current.B19IC_mag)
    legend({'B01IA_mag','B01IB_mag','B01IC_mag','B19IA_mag','B19IB_mag','B19IC_mag'})

    figure
    hold on
    plot(Voltage.Time,Voltage.B24V_mag)
    plot(Voltage.Time,Voltage.B28V_mag)
    plot(Voltage.Time,Voltage.B29V_mag)
    legend({'B24V_mag','B26V_mag','B27V_mag'})


    figure
    hold on
    plot(Voltage.Time,Syncro.UP_V_ph)
    plot(Voltage.Time,Syncro.DOWN_V_ph)
    legend({'UP_V_ph','DOWN_V_ph'})

    figure
    hold on
    plot(Voltage.Time,Syncro.UP_V_mag)
    plot(Voltage.Time,Syncro.DOWN_V_mag)
    legend({'UP_V_mag','DOWN_V_mag'})

    figure
    hold on
    plot(Voltage.Time,Syncro.f_upstream)
    plot(Voltage.Time,Syncro.f_downstream)
    plot(Voltage.Time,50*ones(length(Syncro.f_upstream)))
    legend({'f_upstream','f_downstream'})

%% update LCL filter

figure
hold on
plot(abs(variables.E(19:21,:))')
plot(abs(solution.E(19:21,:))')

figure
plot((imag(variables.S(19:21,:))' - imag(solution.S(19:21,:))')*Grid_para.A_b)

figure
plot((-imag(variables.S(19:21,:)) ./ (abs(variables.E(19:21,:)).^2))')

omega = 2*pi*50;
G_new = (-(imag(variables.S(19:21,:))-imag(solution.S(19:21,:)) )./ (abs(variables.E(19:21,:)).^2))/omega*Grid_para.Y_b;


figure
plot(G_new')