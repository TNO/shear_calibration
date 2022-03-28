% Inverse design: calculate the effect intensity to a given resistance
%
%
%NOTES:
% fzero is an overkill but at least it can accomodate more general formats


function Prob = inv_design(free_par, fix_par, Prob, Prob_actions, Options, load_comb_ii)
%function [Prob,VR] = inv_design(free_par, fix_par, Prob, Options)    

% -------------------------------------------------------------------------
% Initialization
% -------------------------------------------------------------------------
resistance_model    = Options.resistance_model;
load_combination    = Options.load_combination;
consider_VRmin      = Options.consider_VRmin;
K_FI_repr           = Options.K_FI_repr;

chi1                = fix_par(1);
chi2                = fix_par(2);

theta_R             = Prob.theta_R.repr;
f_cc                = Prob.f_cc.repr;
d                   = Prob.d.repr;
b                   = Prob.b.repr;
Asl                 = Prob.Asl.repr;

d_lower             = Prob.d_lower.repr;
a_to_d_ratio        = Prob.d_lower.repr;

gamma_G             = Prob.G.gamma;
ksi                 = Prob.ksi.mean;

switch lower(load_comb_ii)
    case 'traffic'
        Prob.Q1    = Prob_actions.T;
        Prob.K_Q1  = Prob_actions.K_T;
        Prob.psi01 = Prob_actions.psi_T;
    case 'snow_wind'
        Prob.Q1    = Prob_actions.S;
        Prob.K_Q1  = Prob_actions.K_S;
        Prob.psi01 = Prob_actions.psi_S;
        Prob.Q2    = Prob_actions.W;
        Prob.K_Q2  = Prob_actions.K_W;
        Prob.psi02 = Prob_actions.psi_W;
    case 'snow_imposed'
        Prob.Q1    = Prob_actions.S;
        Prob.K_Q1  = Prob_actions.K_S;
        Prob.psi01 = Prob_actions.psi_S;
        Prob.Q2    = Prob_actions.I;
        Prob.K_Q2  = Prob_actions.K_I;
        Prob.psi02 = Prob_actions.psi_I;
    case 'wind_imposed'
        Prob.Q1    = Prob_actions.W;
        Prob.K_Q1  = Prob_actions.K_W;
        Prob.psi01 = Prob_actions.psi_W;
        Prob.Q2    = Prob_actions.I;
        Prob.K_Q2  = Prob_actions.K_I;
        Prob.psi02 = Prob_actions.psi_I;
    otherwise
        error(['Unknown action: ', load_comb_ii])
end

gamma_Q1            = Prob.Q1.gamma;
psi01               = Prob.psi01.mean;
gamma_Q2            = Prob.Q2.gamma;
psi02               = Prob.psi02.mean; 


% -------------------------------------------------------------------------
% Resistance
% -------------------------------------------------------------------------
switch lower(resistance_model)
    case 'ec2_codified_2019'
        gamma_R = free_par(1);
        VR      = EC2_codified_2019(f_cc, Asl, b, d, theta_R, gamma_R, consider_VRmin);
    case 'ec2_pre_2021'
        gamma_R = free_par(1);
        VR      = EC2_pre_2021(f_cc, Asl, b, d, d_lower, a_to_d_ratio, theta_R, gamma_R, consider_VRmin);
    case 'mc2010_level_ii_codified_2019'
        gamma_R = free_par(1);
        VR      = MC2010_level_II_codified_2019(f_cc, Asl, b, d, d_lower, a_to_d_ratio, theta_R, gamma_R);       
    otherwise
        error(['Unknown resistance model:', resistance_model])
end
        
% -------------------------------------------------------------------------
% Effect
% -------------------------------------------------------------------------
switch lower(load_combination)
    case 'ec2_simple'
        % the order of Q1 and Q2 does NOT matter (taken care by the load comb.
        % function)
        VE      = @(x) simple_load_comb(gamma_G, x, gamma_Q1, psi01, x.*chi1/(1-chi1), gamma_Q2, psi02, x.*chi2/(1-chi2));
    case 'ec2_advanced'
        % the order of Q1 and Q2 does NOT matter (taken care by the load comb.
        % function)
        VE      = @(x) advanced_load_comb(gamma_G, ksi, x, gamma_Q1, psi01, x.*chi1/(1-chi1), gamma_Q2, psi02, x.*chi2/(1-chi2));    
    otherwise
        error(['Unknown load combination rule:', load_combination])
end

% -------------------------------------------------------------------------
% Design
% -------------------------------------------------------------------------
% full utilization, E_d = R_d
G_k             = fzero(@(x) VR - K_FI_repr * VE(x), 30);
Q1_k            = G_k*chi1/(1-chi1);
Q2_k            = G_k*chi2/(1-chi2);

% -------------------------------------------------------------------------
% Post-process & collect results
% -------------------------------------------------------------------------
[~, Prob.G]     = Prob.G.repr2mean(G_k, Prob.G);

[~, Prob.Q1]     = Prob.Q1.repr2mean(Q1_k, Prob.Q1);

if ~strcmp(load_comb_ii, 'traffic')
    [~, Prob.Q2]     = Prob.Q2.repr2mean(Q2_k, Prob.Q2);
end

end