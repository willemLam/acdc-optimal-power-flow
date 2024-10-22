%%%%%%%%%%%%%%
%% Get data %%
%%%%%%%%%%%%%%
clear all

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

%% try 1 -> 8 times
% FileName = fullfile(['C:\Users\admin\Documents\Willem\acdc_centralized_controller\results\07-01-2024 15-46_island_EXPERIMENT.mat']);

%% try 3 -> 8 times
% FileName = fullfile(['C:\Users\admin\Documents\Willem\acdc_centralized_controller\results\07-02-2024 09-16_island_EXPERIMENT.mat']);

%% try 3 -> 2afe
% FileName = fullfile(['C:\Users\admin\Documents\Willem\acdc_centralized_controller\results\07-02-2024 13-12_island_EXPERIMENT.mat']);
%start at 6000

%% try 4 -> 2afe GOOD
% FileName = fullfile(['C:\Users\admin\Documents\Willem\acdc_centralized_controller\results\07-03-2024 17-49_island_EXPERIMENT.mat']);


%% try 5 -> 2afe + PV
FileName = fullfile(['C:\Users\admin\Documents\Willem\acdc_centralized_controller\results\07-11-2024 17-07_connected_EXPERIMENT.mat']);



% results\07-Sep-2023_EXPERIMENT_6.mat

% results\08-Sep-2023_EXPERIMENT_6.mat
% results\11-Sep-2023_EXPERIMENT_1.mat
% results\11-Sep-2023_SIMULATION_1.mat

s = whos('-file',FileName);
N = length({s.name});

for i = 1:N
name = ['data_',num2str(i)]
D = load(FileName, genvarname(name));
d= D.(genvarname(name));

    Edc(:,i) = d.Edc_raw;
    Eac(:,i) = d.Eac_raw;
    Iac(:,i) = d.Iac_raw;
    Idc(:,i) = d.Idc_raw;

end
   
%%%%%%%%%%%%%%%%
%% make plots %%
%%%%%%%%%%%%%%%%

f2 = figure(10)
subplot(2,2,1)
hold on
title('Voltage VDC')
plot(Edc','Linewidth',2)
ylabel('Vdc [V]')
xlabel('time')
% xlim([0,3500])

subplot(2,2,2)
hold on
title('Current IDC')
plot(Idc'/10,'Linewidth',2)
ylabel('Idc [A]')
xlabel('time')
% xlim([2600,3500])
hold off

subplot(2,2,3)
hold on
title('Voltage VAC')
plot(abs(Eac(1:3,:))','Linewidth',2)
ylabel('Vac [A]')
xlabel('time')
% xlim([2600,3500])
hold off

subplot(2,2,4)
hold on
title('Current IAC')
plot(abs(Iac(1:3,:))'/10,'Linewidth',2)
ylabel('Iac [A]')
xlabel('time')
% xlim([2600,3500])
hold off

%%

figure
hold on
title('Current IAC')
plot(abs(Iac(:,:))','Linewidth',2)
ylabel('Iac [I]')
xlabel('time')
% xlim([2600,3500])
hold off

figure
hold on
title('Voltage VAC')
plot(abs(Eac(:,:))','Linewidth',2)
ylabel('Vac [V]')
xlabel('time')
% xlim([2600,3500])
hold off

%%

Eac(1)




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

