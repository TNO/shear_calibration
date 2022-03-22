% MC2010 level II shear resistance formula
%
% VR = MC2010_LEVEL_II_CODIFIED_2019(fc, Asl, b, d, C_c, gamma_C)
%
% MIND THE UNITS! The formula is dimensionally inconsistent.
%
%INPUT
% fc        concrete compressive strength, [MPa]
% Asl       area of tensile reinforcement in the considered section, [mm^2]
% b         width, [mm]
% d         depth, [mm]
% C_c       model uncertainty factor, [-]
% gamma_C   partial factor, [-]
%
%
%OUTPUT
% VR    shear resistance, [kN]

function VR = MC2010_level_II_codified_2019(fc, Asl, b, d, C_c, gamma_C)

% -------------------------------------------------------------------------
% MC2010 - level II, according to par. 7.3.3.2
% -------------------------------------------------------------------------
z           = 0.9 .* d;
adratio     = 3.0;                  % treat as input par., later
a           = adratio .* d;
Es          = 2.1e5;    
dg          = 16;                   % treat as input par., later

kdg         = 32 ./ (16 + dg);      % according to eq. (7.3-20) 
kdg(kdg < 0.75) = 0.75;

% calculate design shear resistance in [kN], according to eqs. (7.3-17) and (7.3-21):
V0          = 2.5e3;
n           = length(fc);
VR          = nan(n,1);
for ii = 1:n
    epsx  = @(V) ( (V*a(ii) / z(ii)) + V ) / (2*Es*Asl(ii));
    kv    = @(V) (0.4 / (1 + (1500 * epsx(V)))) * (1300 / (1000 + (kdg*z(ii))));
    VR(ii) = fzero(@(V) (C_c * kv(V) * (sqrt(fc(ii)) / gamma_C) * b(ii) * z(ii)) - V, V0);
end
VR          = 1e-3 .* VR;

% calculate max. allowable total design shear resistance in [kN], according to eq. (7.3-26):
% vmin            = 0.035 .* (k.^1.5) .* (fc.^0.5);   % according to eq. (6.3N)
% VRmin           = 1e-3 .* (vmin .* b .* d);
% idx             = VR - VRmin < 0;
% VR(idx)         = VRmin(idx);
% 
% ID(idx)         = 2;

end