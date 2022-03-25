% MC2010 level II shear resistance formula.
%
% VR = MC2010_LEVEL_II_CODIFIED_2019(fc, Asl, b, d, dg, C_c, gamma_C)
%
% MIND THE UNITS! The formula is dimensionally inconsistent.
% The implementation is based on `Model Code 2010. Final draft. September
% 2011`. We added the `C_c` parameter, it is not part of the MC2010 model
% (`C_c = 1.0` leads to the codified MC2010 model).
%
%INPUT
% fc        concrete compressive strength, [MPa]
% Asl       area of tensile reinforcement in the considered section, [mm^2]
% b         width, [mm]
% d         depth, [mm]
% dg        maximum aggregate size [mm]
% C_c       model uncertainty factor, [-]
% gamma_C   partial factor, [-]
%
%
%OUTPUT
% VR    shear resistance, [kN]

function VR = MC2010_level_II_codified_2019(fc, Asl, b, d, dg, C_c, gamma_C)

% -------------------------------------------------------------------------
% MC2010 - level II, according to par. 7.3.3.2
% -------------------------------------------------------------------------
z           = 0.9 .* d;
adratio     = 3.0;                  % treat as input par., later
a           = adratio .* d;
Es          = 2.1e5;                % [N/mm2]

dg_for_kdg = dg;
% "to account for the loss of aggregate interlock in the cracks due to
% fracture of aggregate particles."
dg_for_kdg(fc > 70) = 0;

kdg         = 32 ./ (16 + dg_for_kdg);   % according to eq. (7.3-20)
kdg(kdg < 0.75) = 0.75;

% calculate design shear resistance in [kN], according to eqs. (7.3-17) and (7.3-21):
V0          = 2.5e3;
n           = length(fc);
VR          = nan(n,1);
for ii = 1:n
    epsx  = @(V) epsx_fun(V, a(ii), z(ii), Asl(ii));
    % eq. (7.3-21)
    kv    = @(V) (0.4 / (1 + (1500 * epsx(V)))) * (1300 / (1000 + (kdg(ii) * z(ii))));
    VR(ii) = fzero(@(V) (C_c * kv(V) * (min(sqrt(fc(ii)), 8) / gamma_C) * b(ii) * z(ii)) - V, V0);
end
VR          = 1e-3 .* VR;


% ..................................
% Utility function
% ..................................

    function e = epsx_fun(V, a, z, Asl)
        % "longitudinal strain is calculated at the mid-depth of the
        % effective shear depth or core layer"
        % based on eq. (7.3-16) and the accompanying text
        e = ((V * a / z) + V ) / (2 * Es * Asl);
        e(e > 0.003) = 0.003;
        e(e < 0) = 0;
    end

end