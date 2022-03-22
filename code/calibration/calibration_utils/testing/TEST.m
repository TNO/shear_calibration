clearvars; close all; clc;


%% INPUT
DS_resist           = 3.5;      % design resistance value coming from design scenario
chi1                = 0.3;
chi2                = 0.0;

beta_target         = 4.7;

ksi                 = 0.85;
psi01               = 0.70;
psi02               = 0.00;

% Resistance [kN]
Prob.R.mean         = 1;
Prob.R.cov          = 0.10;
Prob.R.dist         = 2;
Prob.R.P_char       = NaN;
Prob.R.gamma        = NaN;

% Permanent load [kN]
Prob.G.mean         = 1;
Prob.G.cov          = 0.10;
Prob.G.dist         = 1;
Prob.G.P_char       = 0.5;
Prob.G.gamma        = 1.35;

% Imposed load (annual maxima) - intensity [kN]
Prob.Q1.mean        = 1.0;
Prob.Q1.cov         = 0.53;
Prob.Q1.dist        = 11;
Prob.Q1.P_char      = 0.98;
Prob.Q1.gamma       = 1.50;

Prob               = update_Prob(Prob);


%% INVERSE DESIGN
gamma_G  = Prob.G.gamma;
gamma_Q1 = Prob.Q1.gamma;
gamma_Q2 = 0;

VR       = DS_resist;

load_combination = 'ec2-simple';

switch lower(load_combination)
    case 'ec2-simple'
        % the order of Q1 and Q2 does NOT matter (taken care by the load comb.
        % function)
        VE      = @(x) simple_load_comb(gamma_G, x, gamma_Q1, psi01, x.*chi1/(1-chi1), gamma_Q2, psi02, x.*chi2/(1-chi2));
    case 'ec2-advanced'
        % the order of Q1 and Q2 does NOT matter (taken care by the load comb.
        % function)
        VE      = @(x) advanced_load_comb(gamma_G, ksi, x, gamma_Q1, psi01, x.*chi1/(1-chi1), gamma_Q2, psi02, x.*chi2/(1-chi2));    
    otherwise
        error(['Unknown load combination rule:', load_combination])
end

% full utilization, E_d = R_d
G_k             = fzero(@(x) VR - VE(x), 30);
Q1_k            = G_k*chi1/(1-chi1);
Q2_k            = G_k*chi2/(1-chi2);

Prob.G          = char2mean(G_k, Prob.G);
Prob.Q1         = char2mean(Q1_k, Prob.Q1);
