function [ GNI_estimated, cosIA ] = getTransposedIrradiance( ...
 G, G_cs, B_cs,  panel_slope, panel_azimuth,  ZENITH, AZIMUTH )


incidenceAngle = ...
    (acos(...
    cos(deg2rad(ZENITH)) * cos(deg2rad(panel_slope))  + ...
    sin(deg2rad(panel_slope)) * sin(deg2rad(ZENITH)) .* cos(deg2rad(AZIMUTH - panel_azimuth))));

cosIA = cos(incidenceAngle);
%%%%%%%%%%%
% Decompositon of all-sky irradiance into direct and diffusive components

kc = G./G_cs;
kc = max(kc, 0.1);
% eq. 17
kcmax = 27.21*exp(-114*cos(deg2rad(ZENITH))) + 1.665*exp(-4.494*cos(deg2rad(ZENITH))) + 1.08;
kc = min(kc, kcmax);


% Determine beam real-sky radiation.
B = nan(size(G_cs));
% Eq. 18
B(kc < 16/69) = 0;
B(kc > 1) = B_cs(kc > 1).*kc(kc>1);
sel = (kc >= 16/69 & kc <= 1);
B(sel) = B_cs(sel) .* (kc(sel) -0.38*(1-kc(sel)).^2.5);

% Diffuse real-sky irradiance
D = G - B;

% Transposition to an arbitrary tilted plane
F = 1 - (D./(B + D)).^2;


% global normal irradiance
GNI_estimated = 1 * D .* ((1 + cos(deg2rad(panel_slope)))/2) .* ...
    (1 + F .* sin(deg2rad(panel_slope)/2).^3).* ...
    (1 + F .* cosIA.^2 .* (sin(deg2rad(ZENITH))).^3)  + ...
    1 * B .* (cosIA./cos(deg2rad(ZENITH)));

GNI_estimated = max(GNI_estimated, zeros(size(GNI_estimated)));

end

