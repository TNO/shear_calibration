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

function g = g_fun_Ck(C, f_cc, d, b, Asl, V_Rc_repr, resi_model)

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
% Performance function
% -------------------------------------------------------------------------
g = VR - V_Rc_repr;

end