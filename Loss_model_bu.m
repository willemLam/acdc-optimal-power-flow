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
FilePath = 'C:\Users\admin\Documents\Willem\acdc_centralized_controller\results\';
% FileName_island = fullfile([FilePath, '08-27-2024 16-55_EXPERIMENT_island']);
% FileName_connected = fullfile([FilePath, '08-27-2024 16-55_EXPERIMENT_connected']);
% range = 3:157;
% 
% solution_island = Get_all_data(FileName_island);
% solution_connected = Get_all_data(FileName_connected);
% 
% % Define the modes for each table
% modes_connected = {'Grid connected', 'Prepare Island', 'Forming but connected'};
% modes_island = {'Island', 'Prepare for connection'};
% 
% filtered_connected = solution_connected(ismember(solution_connected.mode, modes_connected), :);
% filtered_island = solution_island(ismember(solution_island.mode, modes_island), :);
% 
% solution = union(filtered_connected,filtered_island);


FileName = fullfile([FilePath, '08-26-2024 17-37_EXPERIMENT_connected']);
FileName = fullfile([FilePath, '08-26-2024 17-37_EXPERIMENT_island']);
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
        solution.feasible(i) = d.proceed;
        solution.duration(i) = d.duration;
        solution.mode{i} = d.mode;
    
        variables.E(:,i) = d.variables.E;
        variables.S(:,i) = d.variables.S;
        variables.I(:,i) = d.variables.I;
    end


%% Visualize data
%power in B19 and B22
figure
hold on
plot(real(variables.S(19,:)))
plot(real(variables.S(22,:)))
plot(real(variables.S(22,:)) + real(variables.S(19,:)))
legend('B19','B22','Delta')

%power in B15
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
Y4 = I_sw./(variables.E(22:24,:)*Grid_para.Vdc_b*Grid_para.Ydc_b);
P_sw_dc = real(variables.E(22:24,:)).*real(variables.E(22:24,:)).*Y4;

I_sw = 0.006*abs(variables.I(19:21,:))*Grid_para.I_b;
Y4 = I_sw./(variables.E(22:24,:)*Grid_para.Vdc_b*Grid_para.Ydc_b);
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
V_diode_forward_dc = interp1(DIODE_piecewise(:,1),DIODE_piecewise(:,2),Imag_dc*Grid_para.Idc_b)./Grid_para.V_b; %pu
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
V_igbt_forward_dc = interp1(IGBT_piecewise(:,1),IGBT_piecewise(:,2),Imag_dc*Grid_para.Idc_b)./Grid_para.V_b; %pu
Zigbt_c_dc = 4/pi*V_igbt_forward_dc./Imag_dc;
V_igbt_forward_ac = interp1(IGBT_piecewise(:,1),IGBT_piecewise(:,2),Imag_ac*Grid_para.I_b)./Grid_para.V_b; %pu
Zigbt_c_ac = 4/pi*V_igbt_forward_ac./Imag_ac;
        
P_c_dc = real((variables.I(22:24,:).*(Zdiode_c_dc+Zigbt_c_dc)) .* conj(variables.I(22:24,:))); %22:24
P_c_ac = real((variables.I(22:24,:).*(Zdiode_c_ac+Zigbt_c_ac)) .* conj(variables.I(22:24,:)));


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
plot((140*P_sw_ac(3,range)' + 10*P_c_dc(3,range)')*Grid_para.A_b)
% plot((9*P_sw_ac(1,range)' + 400*P_c_dc(1,range)')*Grid_para.A_b)
ylabel('loss IC1 W')
legend('observed loss','reconstructed loss')
hold off

subplot(3,1,2)
hold on
plot(abs(real(variables.S(20,range))+real(variables.S(23,range))*Grid_para.A_b))
plot((1*P_sw_ac(1,range)' + 250*P_c_dc(1,range)')*Grid_para.A_b)
% plot((1*P_sw_ac(2,range)' + 180*P_c_dc(2,range)')*Grid_para.A_b)
ylabel('loss IC2 W')
legend('observed loss','reconstructed loss')
hold off

subplot(3,1,3)
hold on
plot(abs(real(variables.S(21,range))+real(variables.S(24,range))*Grid_para.A_b))
plot((150*P_sw_dc(2,range)' + 150* P_c_ac(2,range)')*Grid_para.A_b)
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

%% Voltages
range = 1:390;
figure
subplot(2,1,1)
hold on
plot(abs(abs(variables.E(19,range)))*Grid_para.V_b)
plot(abs(abs(variables.E(20,range)))*Grid_para.V_b)
plot(abs(abs(variables.E(21,range)))*Grid_para.V_b)

hold off

subplot(2,1,2)
hold on
plot(abs(real(variables.E(22,range)))*Grid_para.Vdc_b)
plot(abs(real(variables.E(23,range)))*Grid_para.Vdc_b)
plot(abs(real(variables.E(24,range)))*Grid_para.Vdc_b)
hold off
