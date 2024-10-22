%%%%%%%%%%%%%%
%% Get data %%
%%%%%%%%%%%%%%
clear all
close all


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
FileName = fullfile(['C:\Users\admin\Documents\Willem\acdc_centralized_controller\results\07-02-2024 09-16_island_EXPERIMENT.mat']);
range = 1:700;

%% try 3 -> 2afe
% FileName = fullfile(['C:\Users\admin\Documents\Willem\acdc_centralized_controller\results\07-02-2024 13-12_island_EXPERIMENT.mat']);
% range = 6650:6910

%% try 4 -> 2afe GOOD
% FileName = fullfile(['C:\Users\admin\Documents\Willem\acdc_centralized_controller\results\07-12-2024 14-10_island_EXPERIMENT.mat']);
% range = 100:555;

%% try 5 -> 2afe GOOD
% FileName = fullfile(['C:\Users\admin\Documents\Willem\acdc_centralized_controller\results\07-22-2024 17-28_connected_EXPERIMENT.mat']);

%% try 6 -> 2afe samsung
FileName = fullfile(['C:\Users\admin\Documents\Willem\acdc_centralized_controller\results\07-23-2024 16-44_island_EXPERIMENT.mat']);
range = 120:750;

%% try 7 -> 2afe car
%FileName = fullfile(['C:\Users\admin\Documents\Willem\acdc_centralized_controller\results\07-23-2024 17-55_island_EXPERIMENT.mat']);
% range = 1:825;


%% try 8 -> 2afe 
FileName = fullfile(['C:\Users\admin\Documents\Willem\acdc_centralized_controller\results\07-24-2024 17-18_island_EXPERIMENT.mat']);
% range = 1:825;
% FileName = fullfile(['C:\Users\admin\Documents\Willem\acdc_centralized_controller\results\07-24-2024 17-38_island_EXPERIMENT.mat']);
% range = 1:825;

s = whos('-file',FileName);
N = length({s.name})

for i = 1:N
name = ['data_',num2str(i)]
D = load(FileName, genvarname(name));
d= D.(genvarname(name));

    Edc_full(:,i) = d.Edc_raw;
    Idc_full(:,i) = d.Idc_raw;

end
   
%%
% 
% range = 1:825;
% Edc = movmean(Edc_full(:,range)',50)';
% Idc = movmean(Idc_full(:,range)',50)';

Edc = movmean(Edc_full(:,:)',5)';
Idc = movmean(Idc_full(:,:)',5)';

%%%%%%%%%%%%%%%%
%% make plots %%
%%%%%%%%%%%%%%%%

f2 = figure(10)
subplot(2,1,1)
hold on
title('Voltage VDC')
plot(Edc','Linewidth',2)
ylabel('Vdc [V]')
xlabel('time')
% xlim([0,3500])

subplot(2,1,2)
hold on
title('Current IDC')
plot(Idc'/10,'Linewidth',2)
ylabel('Idc [A]')
xlabel('time')
% xlim([2600,3500])
hold off

%% estimation


RL_1 = -(Edc([1,2],:) - Edc([5,6],:)) ./ (Idc([1,2],:)/10);
RL_2 = -(Edc([5,6],:) - Edc([1,2],:)) ./ (Idc([5,6],:)/10);

RT_1 = (Edc(6,:) - Edc(5,:)) ./ (Idc(2,:)/10);
RT_2 = (Edc(5,:) - Edc(6,:)) ./ (Idc(6,:)/10);

f2 = figure
subplot(2,1,1)
hold on
title('R')
plot(RL_1(1:2,:)','Linewidth',2)
plot(RL_2(1:2,:)','Linewidth',2)
ylabel('R_l')
xlabel('time')
% xlim([0,3500])

subplot(2,1,2)
hold on
title('R')
plot(RT_1(1,:)','Linewidth',2)
plot(RT_2(1,:)','Linewidth',2)
ylabel('R_t')
xlabel('time')
% xlim([2600,3500])
hold off

mean([RL_1,RL_2]')
mean([RT_1,RT_2]')