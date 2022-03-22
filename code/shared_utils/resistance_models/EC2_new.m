% EC2 based shear resistance formula, new, it will be codified
%
% [VR, ID] = EC2_new(fc, Asl, b, d, C_c, gamma_C, consider_vrmin)
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
% ID    indicator of which formula is "active"/governing:
%           1: VR
%           2: VRmin
%
%TODO: the current name EC2_new is bad. It should be changed to something
%more descriptive and constant in time. I do not know enough about the new
%formula to do this.


function [VR, ID] = EC2_new(fc, Asl, b, d, C_c, gamma_M, gamma_C)

% -------------------------------------------------------------------------
% EC2, according to par. 6.2.2
% -------------------------------------------------------------------------
k               = 1 + sqrt(200 ./ d);    
rho_l           = Asl ./ (b .* d);
k(k > 2.0)      = 2.0; 
rho_l(rho_l > 0.02) = 0.02; 

% calculate design shear resistance in [kN], according to eq. (6.2.a):
VR              = 1e-3 .* (C_c/gamma_M .* k .* ((100 .* rho_l .* (fc/gamma_C)).^(1/3)) .* b .* d);
ID              = ones(size(VR));

% tt              = (k .* ((100 .* rho_l).^(1/3)) .* b .* d);

% calculate min. design shear resistance in [kN], according to eq. (6.2.b):
vmin            = 0.035 .* (k.^1.5) .* (fc.^0.5);   % according to eq. (6.3N)
VRmin           = 1e-3 .* (vmin .* b .* d);
idx             = VR - VRmin < 0;
VR(idx)         = VRmin(idx);

ID(idx)         = 2;

end