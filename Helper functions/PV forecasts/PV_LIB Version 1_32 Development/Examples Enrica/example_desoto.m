% load IV curve data for a 36 cell Mitsubishi cSi module
%load 'Desoto_demo.mat'
%load 'data_ecsolar'
clear all
%load 'Desoto_demo.mat'
load('IVecsolar.mat')
% Build structure for constants
Const.E0 = 1000; % W/m2
Const.T0 = 25; % C
Const.k = 1.38066E-23; % J/K
Const.q = 1.60218E-19; % c

% set control variables for parameter estimation
maxiter = 20;
eps1 =NaN; % use default
graphic = false;
n = NaN;  % algorithm will estimate n from data

% Estimate PVsyst model parameters. Structure Model contains the model
% parameters along with array u (true for IV curves where valid parameter
% sets were determined) and vectors Iph, Io, Rsh, Rs which contain
% parameter values for each IV curve.
[Model] = pvl_desoto_parameter_estimation(IVCurves, Specs, Const, maxiter, eps1, graphic, n);

% Calculate 5 parameters for each IV curve using fitted model
dEgdT = -0.0002677;  % use default temperature coefficient for band gap Eg from [1]
[IL, I0, Rs, Rsh, nNsVth] = pvl_calcparams_desoto([IVCurves.Ee],[IVCurves.Tc],Specs.aIsc,Model,Model.Eg_ref, dEgdT);

% Compute IV curves corresponding to each set of 5 parameters
Modeled = pvl_singlediode(IL, I0, Rs, Rsh, nNsVth);

% Compute IV curve at STC
[IL0, Io0, Rs0, Rsh0, nNsVth0] = pvl_calcparams_desoto(1000,25,Specs.aIsc,Model,Model.Eg_ref, dEgdT);
% Compute IV curve points using 5 parameters. Structure STCModeled includes
% fields Imp, Vmp, Isc, Voc, Pmp which will be modeled values at STC
STCModeled = pvl_singlediode(IL0, Io0, Rs0, Rsh0, nNsVth0);

FF0 = (STCModeled.Imp*STCModeled.Vmp)/(STCModeled.Isc*STCModeled.Voc);  % Fill factor at STC

% Calculate temperature coefficients for Vmp, Imp and Voc from simulated IV
% curves at Ee = 1000. Create vector of cell temperatures for simulation to determine temperature
% coefficients
vTc = [15 25 35 45 55 65]';
% Compute sets of parameter values for simulation to determine temperature
% coefficients
[ILtc, Iotc, Rstc, Rshtc, nNsVthtc] = pvl_calcparams_desoto(1000,vTc,Specs.aIsc,Model,Model.Eg_ref, dEgdT);
% Compute IV curves for simulation to determine temperature
% coefficients
Modeledtc = pvl_singlediode(ILtc, Iotc, Rstc, Rshtc, nNsVthtc);

% extract IV curve points and estimate temperature coefficients by
% regression
col1 = ones(size(vTc));
[beta_Vmp]=[vTc-25 col1]\[Modeledtc.Vmp];
betaVmp = beta_Vmp(1);      % temperature coefficient (V/C) for voltage at maximum power
[alpha_Imp]=[vTc-25 col1]\[Modeledtc.Imp];
alphaImp = alpha_Imp(1);    % temperature coefficient (W/C) for current at maximum power

gammaPmp = alphaImp*STCModeled.Vmp+betaVmp*STCModeled.Imp;  % temperature coefficient for maximum power

% Compare calculated and measured IV curve points

Data.Isc = [IVCurves.Isc];
Data.Imp = [IVCurves.Imp];
Data.Voc = [IVCurves.Voc];
Data.Vmp = [IVCurves.Vmp];
Data.E = [IVCurves.Ee];
Data.Isc = Data.Isc(:);
Data.Imp = Data.Imp(:);
Data.Voc = Data.Voc(:);
Data.Vmp = Data.Vmp(:);
Data.E = Data.E(:);

% Plot comparison
figure('Position',[50 50 50+9*96 50+6*96])

Meas = [Data.Isc, Data.Imp, Data.Imp.*Data.Vmp, Data.Voc, Data.Vmp];
Mod = [Modeled.Isc, Modeled.Imp, Modeled.Pmp, Modeled.Voc, Modeled.Vmp];
Title = {'Isc', 'Imp', 'Pmp', 'Voc', 'Vmp'};


for i = [1:5]
    subplot(2,3,i)
    scatter(Meas(:,i), Mod(:,i), 5, 'k', 'filled')
    hold on
    title(Title(i),'fontsize',12,'fontweight','b')
    xlabel('Measured','fontsize',12)
    ylabel('Modeled','fontsize',12)
    axis square
    datamin = min([Meas(:,i);Mod(:,i)]);
    datamax = max([Meas(:,i);Mod(:,i)]);
    datamin = max(datamin - (datamax-datamin)/10,0);
    datamax = datamax + (datamax-datamin)/10;
    plot([datamin, datamax], [datamin, datamax],'r', 'LineWidth',2)
    axis([datamin datamax datamin  datamax])
    box on

end