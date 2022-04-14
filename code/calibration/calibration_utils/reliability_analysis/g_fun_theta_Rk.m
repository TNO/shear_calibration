%Evaluate the performance function
%
%SYNOPSYS
% g = G_FUN(...)
%
%INPUT
% ...
%
%OUTPUT
% g     performance function value

function g = g_fun_theta_Rk(theta_R, f_cc, d, b, Asl, d_lower, a_to_d_ratio, resi_model, V_Rc_repr, consider_VRmin, consider_VRbase)

% =========================================================================
%  PRE-PROCESSING
% =========================================================================
resistance_model = translate_model(resi_model);

% =========================================================================
%  EVALUATION
% =========================================================================
% -------------------------------------------------------------------------
% Resistance
% -------------------------------------------------------------------------
switch lower(resistance_model)
    case 'ec2_codified_2019'
        gamma_R = 1;
        VR      = EC2_codified_2019(f_cc, Asl, b, d, theta_R, gamma_R, consider_VRmin, consider_VRbase);
    case 'ec2_pre_2021'
        gamma_R = 1;
        VR      = EC2_pre_2021(f_cc, Asl, b, d, d_lower, a_to_d_ratio, theta_R, gamma_R, consider_VRmin);
    case 'mc2010_level_ii_codified_2019'
        gamma_R = 1;
        VR      = MC2010_level_II_codified_2019(f_cc, Asl, b, d, d_lower, a_to_d_ratio, theta_R, gamma_R);     
    otherwise
        error(['Unknown resistance model:', resistance_model])
end
        

% -------------------------------------------------------------------------
% Performance function
% -------------------------------------------------------------------------
g = VR - V_Rc_repr;

end