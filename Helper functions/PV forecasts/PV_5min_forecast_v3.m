%PV power forecast computation for 5 minute resolution
function [Pac_perun, Pac_solarmax] =PV_5min_MPPT_forecast_v3(GHI,Tair, MatlabTime, UTC)
%%

DN = MatlabTime;
GNI_estimated = TransposedIrradiance(GHI, DN, UTC);
% GNI_estimated(1) = GNI_current; 
S = double(GNI_estimated);
T = double(Tair);

nm1 = 51;
nm2 = 28;
nm3 = 43 + 43;

[Pdc_perun] =PVModel(S, T, nm1);
[Pdc_solarmax_old] = PVModel(S, T, nm2);
[Pdc_solarmax_new] = PVModel_solarmax(S, T, nm3);

Pdc_solarmax_net = Pdc_solarmax_new+Pdc_solarmax_old;

% [Pdc2] =PVModel(S, T, nm2);
Pac_perun = 0.8*Pdc_perun;% considering efficiency of conveter
Pac_solarmax = 0.85*Pdc_solarmax_net;

% PV_forecast_hour = (-218.84 + 1.0061.*Pdc1 - 9.5767e-6.*Pdc1.^2); % Interpolated function from Encrica

end