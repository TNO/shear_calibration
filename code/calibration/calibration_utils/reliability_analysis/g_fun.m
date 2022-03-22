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

function g = g_fun(C, f_cc, d, b, Asl, G, K_G, ksi, Q1, K_Q1, psi01, Q2, K_Q2, psi02, K_E, resi_model, load_comb)

% =========================================================================
%  PRE-PROCESSING
% =========================================================================
resistance_model = translate_model(resi_model);
load_combination = translate_model(load_comb);

% =========================================================================
%  EVALUATION
% =========================================================================
% -------------------------------------------------------------------------
% Resistance
% -------------------------------------------------------------------------
switch lower(resistance_model)
    case 'ec2_codified_2019'
        gamma_C = 1;
        VR      = EC2_codified_2019(f_cc, Asl, b, d, C, gamma_C);
    case 'ec2_new'
        gamma_M = 1;
        gamma_C = 1;
        VR      = EC2_new(f_cc, Asl, b, d, C, gamma_M, gamma_C);
    case 'ec2_proposed_tg4_2016'
        gamma_C = 1;
        VR      = EC2_proposed_TG4_2016(f_cc, Asl, b, d, C, gamma_C);
    case 'ec2_proposed_yuguang_2019'
        gamma_C = 1;
        VR      = EC2_proposed_Yuguang_2019(f_cc, Asl, b, d, C, gamma_C);
    case 'mc2010_level_ii_codified_2019'
        gamma_C = 1;
        VR      = MC2010_level_II_codified_2019(f_cc, Asl, b, d, C, gamma_C);  
    case 'mc2010_new'
        % to be fixed later
%         gamma_M = 1;
%         gamma_C = 1;
%         VR      = MC2010levelII_new(f_cc, Asl, b, d, C, gamma_M, gamma_C);         
    otherwise
        error(['Unknown resistance model:', resistance_model])
end
        
% -------------------------------------------------------------------------
% Effect
% -------------------------------------------------------------------------
switch lower(load_combination)
    case 'ec2_simple'
        gamma_G = 1;
        gamma_Q = 1;
        % the order of Q1 and Q2 does NOT matter (taken care by the load comb.
        % function)
        [VE, ~, ~]      = simple_load_comb(gamma_G, K_G.*G, gamma_Q, psi01, K_Q1.*Q1, gamma_Q, psi02, K_Q2.*Q2);
%         [~, VE1, ~]      = simple_load_comb(gamma_G, K_G.*G, gamma_Q, psi01, K_Q1.*Q1, gamma_Q, psi02, K_Q2.*Q2);
%         VE = VE1;
    case 'ec2_advanced'
        gamma_G = 1;
        gamma_Q = 1;
        % the order of Q1 and Q2 does NOT matter (taken care by the load comb.
        % function)
        VE      = advanced_load_comb(gamma_G, ksi, K_G.*G, gamma_Q, psi01, K_Q1.*Q1, gamma_Q, psi02, K_Q2.*Q2); 
    otherwise
        error(['Unknown load combination rule:', load_combination])
end

% -------------------------------------------------------------------------
% Performance function
% -------------------------------------------------------------------------
g = VR - K_E .* VE;

end