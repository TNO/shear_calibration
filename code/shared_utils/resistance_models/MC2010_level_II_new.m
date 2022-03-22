% MC2010 - level II based resistance formula
%
%INPUT
% fc
%
%
%OUTPUT
% VR shear resistance
% ID indicator of which formula is active:
%           1 VR
%           2 VRmax

function [VR, ID, tt] = MC2010levelII_new(fc, Asl, b, d, C_c, gamma_M, gamma_C)

% -------------------------------------------------------------------------
% MC2010 - level II, according to par. 7.3.3.2
% -------------------------------------------------------------------------
z           = 0.9 .* d;
adratio     = 3.0;          % treat as input par., later
a    0       = adratio .* d;
Es          = 2.1e5;    
dg          = 16;           % treat as input par., later

kdg             = 32 ./ (16 + dg);      % according to eq. (7.3-20) 
kdg(kdg < 0.75) = 0.75;

% calculate design shear resistance in [kN], according to eqs. (7.3-17) and (7.3-21):
V0   = 2.5e3;
for i = 1:length(fc)
    epsx  = @(V) ( (V*a(i) / z(i)) + V ) / (2*Es*Asl(i));
    kv    = @(V) (0.4 / (1 + (1500 * epsx(V)))) * (1300 / (1000 + (kdg*z(i))));
    VR(i) = fzero(@(V) ((C_c(i)/gamma_M) * kv(V) * (sqrt(fc(i)) / gamma_C) * b(i) * z(i)) - V, V0);
end
VR   = 1e-3 .* VR;
ID   = ones(size(VR));

tt   = [];

% calculate max. allowable total design shear resistance in [kN], according to eq. (7.3-26):
% vmin            = 0.035 .* (k.^1.5) .* (fc.^0.5);   % according to eq. (6.3N)
% VRmin           = 1e-3 .* (vmin .* b .* d);
% idx             = VR - VRmin < 0;
% VR(idx)         = VRmin(idx);
% 
% ID(idx)         = 2;

end