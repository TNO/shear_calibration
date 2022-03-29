% Wrapper for the considered shear resistance models (formulas).
function [VR, ID] = shear_formula(fc, Asl, b, d, d_lower, a_to_d_ratio, fsy, gamma_S, theta_R, gamma_R, consider_VRmin, model)
    switch lower(model)
        case 'en1992-1-1'
            [VR, ID] = EC2_codified_2019(fc, Asl, b, d, theta_R, gamma_R, consider_VRmin);
        case 'pren1992-1-1'
            [VR, ID] = EC2_pre_2021(fc, Asl, b, d, d_lower, a_to_d_ratio, fsy, gamma_S, theta_R, gamma_R, consider_VRmin);
        case 'mc2010'
            [VR, ID] = MC2010_level_II_codified_2019(fc, Asl, b, d, d_lower, a_to_d_ratio, theta_R, gamma_R);
        otherwise
            error('Unknown model.')
    end
end