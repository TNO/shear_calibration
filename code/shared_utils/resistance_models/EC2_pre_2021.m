% prEC2 based shear resistance formula (2021-October-14)
%
% [VR, ID]  = EC2_pre_2021(fc, Asl, b, d, d_lower, a_to_d_ratio, theta_R, gamma_R, consider_vrmin)
%
% MIND THE UNITS! The formula is dimensionally inconsistent.
% The implementation is based on `prEN1992-1-1 2021 officiële versies
% van de Formal Enquiries. 2021-10-14`.
% 
% Slightly modified (generalized) to make it reasonable from a structural
% reliability point of view. To get back to codified formula for persistent
% and transient deisign situations (Table 4.3):
% * characteristic: gamma_S = 1.0, theta_R = 1.0, gamma_R = 1.0, consider_vrmin = true
% * design: gamma_S = 1.15, theta_R = 1.0, gamma_R = 1.4, consider_vrmin = true
%
%INPUT
% fc        concrete compressive strength, [MPa]
% Asl       area of tensile reinforcement in the considered section, [mm^2]
% b         width, [mm]
% d         depth, [mm]
% d_lower   size of the maximum sieve gird used for aggregates [mm]
% a_to_d_ratio []
% fsy       yield stress of the rebar designed as flexural reinforcment
%           [MPa]
% gamma_S   rebar material strength partial factor
% theta_R   model uncertainty factor, [-]
% gamma_R   model uncertainty partial factor, [-]
%
%OPTIONAL
% consider_VRmin    Consider VRmin in the calculation?, default: no
%
%OUTPUT
% VR    shear resistance, [kN]
% ID    indicator of which formula is "active"/governing (for this 
%       formula the ID is always 1, ID is included for consistency with
%       other shear formulas):
%           1: VR 
%           2: VRmin


function [VR, ID] = EC2_pre_2021(fc, Asl, b, d, d_lower, a_to_d_ratio, fsy, gamma_S, theta_R, gamma_R, consider_VRmin)

% -------------------------------------------------------------------------
% Initialize
% -------------------------------------------------------------------------
% TODO: make it more general and robust

if length(d_lower) == 1
    d_lower = d_lower * ones(size(fc));
end
if length(a_to_d_ratio) == 1
    a_to_d_ratio = a_to_d_ratio * ones(size(fc));
end
if length(fsy) == 1
    fsy = fsy * ones(size(fc));
end
if length(theta_R) == 1
    theta_R = theta_R * ones(size(fc));
end

% -------------------------------------------------------------------------
% prEN1992-1-1
% -------------------------------------------------------------------------
% lever arm
z               = 0.9 * d;

% size parameter describing the failure zone roughness
% section 8.2.1, paragraph (4)
d_dg            = 16 + d_lower;

idx             = fc > 60;
d_dg(idx)       = 16 + d_lower(idx) .* (60 ./ fc(idx)).^4;
d_dg            = bsxfun(@min, d_dg, 40);

rho_l           = Asl ./ (b .* d); 

% section 8.2.2, paragraph (3)
% effective shear span
a_cs            = a_to_d_ratio .* d;
d_for_tau_Rc_base = d;
idx = d < 4 * a_cs;
d_for_tau_Rc_base(idx) = max(a_cs(idx), d(idx));

% eq. (8.16):
tau_Rc_base     = 0.66 / gamma_R .* ((100 .* rho_l .* fc .* d_dg./d_for_tau_Rc_base).^(1/3));

tau_Rc          = tau_Rc_base;
ID              = ones(size(tau_Rc));

if consider_VRmin == 1
    % Eq.(8.11)
    tau_Rc_min      = 11 / gamma_R * sqrt(fc ./ (fsy ./ gamma_S) .* d_dg ./ d);
    idx             = tau_Rc_base < tau_Rc_min;
    tau_Rc(idx)     = tau_Rc_min(idx);
    ID(idx)         = 2;
end

% Eq.(8.10a)
% 1e-3: [N] -> [kN]
VR = 1e-3 .* theta_R .* tau_Rc .* b .* z;

end