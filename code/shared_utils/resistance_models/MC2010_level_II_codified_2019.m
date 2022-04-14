% MC2010 level II shear resistance formula.
%
% [VR, ID] = MC2010_LEVEL_II_CODIFIED_2019(fc, Asl, b, d, d_lower, a_to_d_ratio, theta_R, gamma_R)
%
% MIND THE UNITS! The formula is dimensionally inconsistent.
% The implementation is based on `Model Code 2010. Final draft. September
% 2011`. We added the `theta_R` parameter, it is not part of the MC2010 model
% (`theta_R = 1.0` leads to the codified MC2010 model).
%
%INPUT
% fc        concrete compressive strength, [MPa]
% Asl       area of tensile reinforcement in the considered section, [mm^2]
% b         width, [mm]
% d         depth, [mm]
% d_lower   size of the maximum sieve gird used for aggregates [mm]
% a_to_d_ratio []
% theta_R   model uncertainty factor, [-]
% gamma_R   partial factor, [-]
%
%
%OUTPUT
% VR    shear resistance, [kN]
% ID    indicator of which formula is "active"/governing (for this 
%       formula the ID is always 1, ID is included for consistency with
%       other shear formulas):
%           1: VR 
%           2: VRmin

function [VR, ID] = MC2010_level_II_codified_2019(fc, Asl, b, d, d_lower, a_to_d_ratio, theta_R, gamma_R)

% -------------------------------------------------------------------------
% Initialize
% -------------------------------------------------------------------------
[fc, Asl, b, d, d_lower, a_to_d_ratio, theta_R, gamma_R] = ...
    equalize_vector_lengths( ...
        fc, Asl, b, d, d_lower, a_to_d_ratio, theta_R, gamma_R ...
    );

% -------------------------------------------------------------------------
% MC2010 - level II, according to par. 7.3.3.2
% -------------------------------------------------------------------------
z = 0.9 .* d;
a = a_to_d_ratio .* d;
Es = 2.1e5;                % [N/mm2]

% dg - maximum aggregate size [mm]
% according to Yuguang's email of 2022-Mar-27
dg = d_lower;

dg_for_kdg = dg;
% "to account for the loss of aggregate interlock in the cracks due to
% fracture of aggregate particles."
dg_for_kdg(fc > 70) = 0;

kdg = 32 ./ (16 + dg_for_kdg);   % according to eq. (7.3-20)
kdg(kdg < 0.75) = 0.75;

% calculate design shear resistance in [kN], according to eqs. (7.3-17) and 

% auxiliary variables
% V independent term in eq. (7.3-21)
c1 = 0.4 * 1300 ./ (1000 + kdg .* z);

% eps_x = V * c2; based on eq. (7.3-16) and the accompanying text
c2 = (a ./ z + 1) ./ (2 * Es .* Asl);

% right hand side of eq. (7.3-17) without kv
c3 = theta_R .* min(sqrt(fc), 8) ./ gamma_R .* b .* z;

% eq. (7.3-17) with the c coefficients: V = c1 * c3 / (1 + 1500*c3*V)
% coeffs of a quadratic equation of the form of: a*x^2 + b*x + c
a = 1500 * c2;
b = 1;
c = -c1 .* c3;

VR = (-b + sqrt(b.^2 - 4*a.*c))./(2*a);

% correct if eps_x limits are reached
eps_x = VR .* c2;
% eps_x < 0; impossible under our assumptions (c2) so it is commented out
% bm_negative_epsx = eps_x < 0;
% VR(bm_negative_epsx) = c1(bm_negative_epsx) .* c3(bm_negative_epsx);

% eps_x > 0.003
bm_large_epsx = eps_x > 0.003;
% 5.5 = 1 + 1500 * 0.003
if any(bm_large_epsx)
    VR(bm_large_epsx) = c1(bm_large_epsx) / 5.5 .* c3(bm_large_epsx);
end

VR = 1e-3 .* VR;

ID = ones(size(VR));

end
