%PV power forecast computation for 5 minute resolution
function [Pac_perun, Pac_solarmax, Pac_facade] =PV_production(GHI,GI, GVI, Tair, MatlabTime, UTC)
%%

DN = MatlabTime;

panel_slope = 10;
panel_slope_facade = 90;
GNI_estimated = GI;%TransposedIrradiance(GHI, DN, UTC,panel_slope);
GNI_estimated_facade = GVI;%TransposedIrradiance(GHI, DN, UTC,panel_slope_facade);
% GNI_estimated(1) = GNI_current; 
S = double(GNI_estimated);
S_facade = double(GNI_estimated_facade);
T = double(Tair);

nm1 = 51;
nm2 = 28;
nm3 = 43;
nm4 = 76;
[Pdc_perun] =PVModel_perun(S, T, nm1);
[Pdc_solarmax_old] = PVModel_solarmax_old(S, T, nm2);
[Pdc_solarmax_new] = PVModel_solarmax_new(S, T, nm3);
[Pdc_facade] = PVModel_facade(S_facade, T, nm4);


Pdc_solarmax_net = Pdc_solarmax_new+Pdc_solarmax_old;

% [Pdc2] =PVModel(S, T, nm2);
Pac_perun = 0.8*Pdc_perun;% considering efficiency of conveter
Pac_solarmax = 0.95*Pdc_solarmax_net;
Pac_facade = 0.95*Pdc_facade;

% PV_forecast_hour = (-218.84 + 1.0061.*Pdc1 - 9.5767e-6.*Pdc1.^2); % Interpolated function from Encrica

end