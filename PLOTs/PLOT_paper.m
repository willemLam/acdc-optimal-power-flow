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

%get the data
% FilePath = 'C:\Users\admin\Documents\Willem\acdc_centralized_controller\results\';
% FileName_island = fullfile([FilePath, '08-27-2024 16-55_EXPERIMENT_island']);
% FileName_connected = fullfile([FilePath, '08-27-2024 16-55_EXPERIMENT_connected']);
% range = 3:157;

% good test - no island
% FileName_island = fullfile([FilePath, '09-06-2024 10-12_EXPERIMENT_connected']);
% FileName_connected = fullfile([FilePath, '09-06-2024 10-12_EXPERIMENT_connected']);
% range = 3:240;

% FilePath = 'C:\Users\admin\Documents\Willem\acdc_centralized_controller\results\';
% FileName_island = fullfile([FilePath, '09-24-2024 15-22_EXPERIMENT_island']);
% FileName_connected = fullfile([FilePath, '09-24-2024 15-21_EXPERIMENT_connected']);
% range = 3:405;

% Nb 1 - good - no pv
% FilePath = 'C:\Users\admin\Documents\Willem\acdc_centralized_controller\results\';
% FileName_island = fullfile([FilePath, '09-24-2024 16-15_EXPERIMENT_island']);
% FileName_connected = fullfile([FilePath, '09-24-2024 16-15_EXPERIMENT_connected']);
% range = 5:520;

% Nb 2 - good - low level pv
FilePath = 'C:\Users\admin\Documents\Willem\acdc_centralized_controller\results\';
FileName_island = fullfile([FilePath, '09-24-2024 18-08_EXPERIMENT_island']);
FileName_connected = fullfile([FilePath, '09-24-2024 18-08_EXPERIMENT_connected']);
range = 3:490;



% Nb trip at resync
FilePath = 'C:\Users\admin\Documents\Willem\acdc_centralized_controller\results\';
FileName_island = fullfile([FilePath, '10-02-2024 15-55_EXPERIMENT_island']);
FileName_connected = fullfile([FilePath, '10-02-2024 15-55_EXPERIMENT_connected']);

FileName_island = fullfile([FilePath, '10-02-2024 16-56_EXPERIMENT_island']);
FileName_connected = fullfile([FilePath, '10-02-2024 16-56_EXPERIMENT_connected']);


FileName_island = fullfile([FilePath, '10-02-2024 15-37_EXPERIMENT_island']);
FileName_connected = fullfile([FilePath, '10-02-2024 15-37_EXPERIMENT_connected']);


%good -5kw big overvoltage
FileName_island = fullfile([FilePath, '10-10-2024 15-27_EXPERIMENT_island']);
FileName_connected = fullfile([FilePath, '10-10-2024 15-27_EXPERIMENT_connected']);

FileName_island = fullfile([FilePath, '10-10-2024 15-56_EXPERIMENT_island']);
FileName_connected = fullfile([FilePath, '10-10-2024 15-56_EXPERIMENT_connected']);

% 5x island - resync
FileName_island = fullfile([FilePath, '10-10-2024 17-26_EXPERIMENT_island']);
FileName_connected = fullfile([FilePath, '10-10-2024 17-26_EXPERIMENT_connected']);

FileName_island = fullfile([FilePath, '10-11-2024 17-09_EXPERIMENT_island']);
FileName_connected = fullfile([FilePath, '10-11-2024 17-09_EXPERIMENT_connected']);

% 10-15: 15:10, 15:27 15:42(time)
FileName_island = fullfile([FilePath, '10-15-2024 15-27_EXPERIMENT_island']);
FileName_connected = fullfile([FilePath, '10-15-2024 15-27_EXPERIMENT_connected']);

% 11:48 -  
FileName_island = fullfile([FilePath, '10-17-2024 12-22_EXPERIMENT_island']);
FileName_connected = fullfile([FilePath, '10-17-2024 12-22_EXPERIMENT_connected']);

% 11:28 -  14:22 - 15:45! - 16:09!
FileName_island = fullfile([FilePath, '10-18-2024 16-09_EXPERIMENT_island']);
FileName_connected = fullfile([FilePath, '10-18-2024 16-09_EXPERIMENT_connected']);

% 11:08 - 11:19 - 11:42 - 11:54 - 13:32! - 13:40! - 13:51! - 13:58! -
% 14:15! 14:28!! 14:37!! 14:53!!! 15:05!!! 15:23!! 15:38!! 15:49 16:07!
FileName_island = fullfile([FilePath, '10-21-2024 15-49_EXPERIMENT_island']);
FileName_connected = fullfile([FilePath, '10-21-2024 15-49_EXPERIMENT_connected']);


% %crash - high level pv
% FilePath = 'C:\Users\admin\Documents\Willem\acdc_centralized_controller\results\';
% FileName_island = fullfile([FilePath, '09-24-2024 18-18_EXPERIMENT_island']);
% FileName_connected = fullfile([FilePath, '09-24-2024 18-18_EXPERIMENT_connected']);
% range = 3:170;

solution_island = Get_all_data(FileName_island);
solution_connected = Get_all_data(FileName_connected);

% Define the modes for each table
modes_connected = {'Grid connected', 'Prepare Island', 'Forming but connected'};
modes_island = {'Island', 'Prepare for connection'};




filtered_connected = solution_connected(ismember(solution_connected.mode, modes_connected), :);
filtered_island = solution_island(ismember(solution_island.mode, modes_island), :);

solution = union(filtered_connected,filtered_island);

% range = 5:min(height(solution_connected),height(solution_island));
range = 5:height(solution);
% color mape for states
modeColors = containers.Map({'Grid connected', 'Prepare Island', 'Forming but connected', 'Island', 'Prepare for connection'}, ...
                                {m_cmap1(1,:),m_cmap3(1,:),m_cmap3(5,:),m_cmap1(5,:),m_cmap3(3,:)});
%%
SAVE = true;

%% PLOT CPU time


f9 = make_plot_time2(solution,solution_connected,solution_island,Grid_para,constraints,range,modeColors);


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