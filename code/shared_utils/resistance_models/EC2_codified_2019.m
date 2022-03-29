% EC2 based shear resistance formula
%
% [VR, ID] = EC2_codified_2019(fc, Asl, b, d, theta_R, gamma_R, consider_vrmin)
%
% MIND THE UNITS! The formula is dimensionally inconsistent.
% Slightly modified (generalized) to make it reasonable from a structural
% reliability point of view. To get back to codified formula:
% * characteristic: theta_R = 1.0, gamma_R = 1.0, consider_vrmin = true
% * design: theta_R = 1.0, gamma_R = 1.5, consider_vrmin = true
%
%INPUT
% fc        concrete compressive strength, [MPa]
% Asl       area of tensile reinforcement in the considered section, [mm^2]
% b         width, [mm]
% d         depth, [mm]
% theta_R   model uncertainty factor, [-]
% gamma_R   partial factor, [-]
%
%OPTIONAL
% consider_VRmin    Consider VRmin in the calculation?, default: no
%
%OUTPUT
% VR    shear resistance, [kN]
% ID    indicator of which formula is "active"/governing:
%           1: VR
%           2: VRmin

function [VR, ID] = EC2_codified_2019(fc, Asl, b, d, theta_R, gamma_R, consider_VRmin)

if nargin < 7
   consider_VRmin = false;
end

% -------------------------------------------------------------------------
% EC2, according to par. 6.2.2
% -------------------------------------------------------------------------
k               = 1 + sqrt(200 ./ d);    
rho_l           = Asl ./ (b .* d);
k(k > 2.0)      = 2.0; 
rho_l(rho_l > 0.02) = 0.02; 

% calculate design shear resistance in [kN], according to eq. (6.2.a):
VR              = 1e-3 .* (0.18 / gamma_R .* k .* ((100 .* rho_l .* fc).^(1/3)) .* b .* d);
ID              = ones(size(VR));

if consider_VRmin == 1
    % calculate min. design shear resistance in [kN], according to eq. (6.2.b):
    % we replaced 0.035 with 0.035 * 1.5 / gamma_C; to make the formula
    % reasonable from a structural reliability point of view
    vmin            = (0.035 * 1.5) / gamma_R .* (k.^1.5) .* (fc.^0.5);   % according to eq. (6.3N)
    VRmin           = 1e-3 .* (vmin .* b .* d);
    idx             = VR < VRmin;
    VR(idx)         = VRmin(idx);
    ID(idx)         = 2;
end

VR = VR .* theta_R;

end