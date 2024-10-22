close all
clear all
addpath(genpath(pwd))

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

% 11-08     : 10A, deltaI=[1.1: 0.5], 3 turn, deltaV=4V, extra:Q not stable, 
% 11-19     : 10A, deltaI=[0.2:-0.1], 1 turn, deltaV=4V,
% 11-42     : 20A, deltaI=[0.8:-0.6], 1 turn, deltaV=8V,
% 11-54     : 20A, deltaI=[1.1:-0.1], 5 turn, deltaV=6V,           extra: sexy, but little spike when preparing for island
% 13-32!    : 20A, deltaI=[0.8:-0.05], 1 turn, deltaV=5V, GCP=10kW, extra: Q not very stable, 
% 13-40!    : 20A, deltaI=[0.8:-0.1], 0 turn, deltaV=8V, GCP=8kW, extra: 
% 13-51!    : 17A, deltaI=[0.8:-0.6], 3 turn, deltaV=6V, GCP=8kW, extra: cloudy
% 13-58!    : 17A, deltaI=[0.9:-0], 5 turn, deltaV=8V, GCP=6kW, extra: Q not very stable, 
% 14-15!    : 17A, deltaI=[1.0:-0.1], 4 turn, deltaV=10V, GCP=8kW, extra: sexy, GCP goes to 0
%-> 14-28!!   : 17A, deltaI=[0.5:-0], 1 turn, deltaV=6V, GCP=8kW, extra: YES, GCP goes to 0
%-> 14-37!!   : 17A, deltaI=[0.5:-0], 1 turn, deltaV=5V, GCP=8kW, extra: YES, GCP goes to 0, little spike when switch opens
%-> 14-53!!!  : 17A, deltaI=[0.3:-0.05], 6 turn, deltaV=8V, GCP=1kW, extra: car, 
% 15-05!!!  : 15A, deltaI=[0.4:-0.1], 2 turn, deltaV=10V, GCP=1kW, extra: car, not very pretty
% 15-23!!   : 15A, deltaI=[0.3:0], 2 turn, deltaV=5V, GCP=0kW, extra:
% 15-38!!   : 15A, deltaI=[0.1:-0.05], 1 turn, deltaV=5V, GCP=6kW, extra:weird resync
% 15-49     : 10A, deltaI=[0.:-0.1],  1 turn, deltaV=2V, GCP=6kW, extra: NOPE! -20Vdc drop
% 16-07!    : NOPE

tst = '14-53'; %'14-28', '14-37', '14-53'
FileName_island = fullfile(['Results/', '10-21-2024 ',tst,'_EXPERIMENT_island']);
FileName_connected = fullfile(['Results/', '10-21-2024 ',tst,'_EXPERIMENT_connected']);



solution_island = Get_all_data(FileName_island);
solution_connected = Get_all_data(FileName_connected);

% Define the modes for each table
modes_connected = {'Grid connected', 'Prepare Island', 'Forming but connected'};
modes_island = {'Island', 'Prepare for connection'};




filtered_connected = solution_connected(ismember(solution_connected.mode, modes_connected), :);
filtered_island = solution_island(ismember(solution_island.mode, modes_island), :);
solution = union(filtered_connected,filtered_island);

range = 5:height(solution);

% color mape for states
modeColors = containers.Map({'Grid connected', 'Prepare Island', 'Forming but connected', 'Island', 'Prepare for connection'}, ...
                                {m_cmap1(1,:),m_cmap3(1,:),m_cmap3(5,:),m_cmap1(5,:),m_cmap3(3,:)});
%
SAVE = false;

%% PLOT CPU time
% f9 = make_plot_time2(solution,solution_connected,solution_island,Grid_para,constraints,range,modeColors);


f0 = make_plot_cputime2(solution,range);
% f0 = make_plot_cputime(solution_island,solution_connected,range);

if SAVE
    folder = './PLOTs/Figures';
    saveas(f0,[folder filesep() 'jpg' filesep() 'Fig-CPU_time'],'jpg');
    saveas(f0,[folder filesep() 'eps' filesep() 'Fig-CPU_time'],'epsc');
end
%% PLOT Reactive power error
f1 = make_plot_Reactive_state_opf(solution,Grid_para,range,modeColors);

if SAVE
    folder = './PLOTs/Figures';
    saveas(f1,[folder filesep() 'jpg' filesep() 'Fig-Reactive_Power_OPF_STATE'],'jpg');
    saveas(f1,[folder filesep() 'eps' filesep() 'Fig-Reactive_Power_OPF_STATE'],'epsc');
end
%% PLOT Active power error
f2 = make_plot_Active_state_opf(solution,Grid_para,range,modeColors);

if SAVE
    folder = './PLOTs/Figures';
    saveas(f2,[folder filesep() 'jpg' filesep() 'Fig-Active_Power_OPF_STATE'],'jpg');
    saveas(f2,[folder filesep() 'eps' filesep() 'Fig-Active_Power_OPF_STATE'],'epsc');
end
%% State of Charge
% f3 = make_plot_SOC(solution,Grid_para,range,modeColors);
f3 = make_plot_SOC_P(solution,Grid_para,range,modeColors);
if SAVE
    folder = './PLOTs/Figures';
    saveas(f3,[folder filesep() 'jpg' filesep() 'Fig-SOC'],'jpg');
    saveas(f3,[folder filesep() 'eps' filesep() 'Fig-SOC'],'epsc');
end
%% Current flow
f4 = make_plot_current(solution,Grid_para,constraints,range,modeColors);

if SAVE
    folder = './PLOTs/Figures';
    saveas(f4,[folder filesep() 'jpg' filesep() 'Fig-Current'],'jpg');
    saveas(f4,[folder filesep() 'eps' filesep() 'Fig-Current'],'epsc');
end

%% DC Voltage
f5 = make_plot_DCVoltage(solution,Grid_para,constraints,range,modeColors);

if SAVE
    folder = './PLOTs/Figures';
    saveas(f5,[folder filesep() 'jpg' filesep() 'Fig-DCVoltage'],'jpg');
    saveas(f5,[folder filesep() 'eps' filesep() 'Fig-DCVoltage'],'epsc');
end

%% AC Voltage
f6 = make_plot_ACVoltage(solution,Grid_para,constraints,range,modeColors);

if SAVE
    folder = './PLOTs/Figures';
    saveas(f6,[folder filesep() 'jpg' filesep() 'Fig-ACVoltage'],'jpg');
    saveas(f6,[folder filesep() 'eps' filesep() 'Fig-ACVoltage'],'epsc');
end

%% freqency
f7 = make_plot_frequency(solution,Grid_para,constraints,range,modeColors);

if SAVE
    folder = './PLOTs/Figures';
    saveas(f7,[folder filesep() 'jpg' filesep() 'Fig-frequency'],'jpg');
    saveas(f7,[folder filesep() 'eps' filesep() 'Fig-frequency'],'epsc');
end

%% angle
f8 = make_plot_angle(solution,Grid_para,constraints,range,modeColors);

if SAVE
    folder = './PLOTs/Figures';
    saveas(f8,[folder filesep() 'jpg' filesep() 'Fig-angle'],'jpg');
    saveas(f8,[folder filesep() 'eps' filesep() 'Fig-angle'],'epsc');
end

%% powerGCP
f9 = make_plot_GCP(solution,Grid_para,constraints,range,modeColors);

if SAVE
    folder = './PLOTs/Figures';
    saveas(f9,[folder filesep() 'jpg' filesep() 'Fig-GCP'],'jpg');
    saveas(f9,[folder filesep() 'eps' filesep() 'Fig-GCP'],'epsc');
end