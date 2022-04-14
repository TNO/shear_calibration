% Get the input data/parameters
%
%SYNOPSYS
% [Prob, Prob_actions, DS] = GET_INPUT(Options)
%
%INPUT
%
%OUTPUT
% Prob          probabilistic models
% Prob_actions  probabilistic models for actions (temp.)
% DS            design scenarios

%NOTES:
% * If you want to consider only a single variable action set the other
%   like this:
%       - zero mean
%       - 0 for distribution type (deterministic/constant)
% * The random variables should be consistent with those of in run_reli.m and
%   related m-files.

%TODO
% * weights!

function [Prob, Prob_actions, DS] = get_input(Options)


% =========================================================================
% DESIGN SCENARIOS
% =========================================================================
resistance_model    = Options.resistance_model;
consider_VRmin      = Options.consider_VRmin;

load_combs          = Options.load_combs;
load_comb_weights   = Options.load_comb_weights;

% effective depth - nominal values
% d_ds                = [150, 300, 450, 600, 750];
% d_ds                = [300, 450];
d_ds                = 300;

% concrete - characteristic compressive strength
% f_cck_ds            = [20, 40, 60, 80];
% f_cck_ds            = [40, 60];
f_cck_ds            = 40;

% load ratio for variable load 1
% chi1_ds         = 0.1:0.01:0.9;
chi1_ds         = 0.1:0.1:0.9;

% load ratio for variable load 2
% chi2_ds         = 0.1:0.2:0.9;
% chi2_ds         = 0.1:0.05:0.9;
chi2_ds         = 0.1:0.1:0.9;

% reinforcement ratio - with nominal values
% rho_ds              = [0.005, 0.010, 0.015];
% rho_ds              = [0.005, 0.010];
rho_ds          = 0.01;

% size of the maximum sieve gird used for aggregates
% d_lower_ds          = [8, 16, 24];
d_lower_ds          = [8, 16];
% d_lower_ds          = 16;

% a-d ratio
% a_to_d_ratio_ds     = [2, 3, 4];
% a_to_d_ratio_ds     = [2, 3];
a_to_d_ratio_ds     = 3;

% ...............................................
% FOR TESTING
% d_ds            = 120;
% f_cck_ds        = 30;
% chi1_ds         = 0.1:0.1:0.9;
% chi2_ds         = 0.3;
% rho_ds          = 0.01;
% ...............................................

DS.Range.d          = d_ds;
DS.Range.f_cck      = f_cck_ds;
DS.Range.chi1       = chi1_ds;
DS.Range.chi2       = chi2_ds;
DS.Range.rho        = rho_ds;
DS.Range.d_lower    = d_lower_ds;
DS.Range.a_to_d_ratio = a_to_d_ratio_ds;

[Chi1q, Chi2q]      = meshgrid(chi1_ds, chi2_ds);

% weights to chi1 and chi2
n_lc                = length(load_combs);
for ii = 1:n_lc
    lc_ii       = load_combs{ii};
    lcw_ii      = load_comb_weights(ii);
    Chi1q_ii    = Chi1q;
    Chi2q_ii    = Chi2q;

    W       = xlsread(Options.weights_filepath, lc_ii, 'B3:Z100', 'basic');
    chi1    = xlsread(Options.weights_filepath, lc_ii, 'A3:A100', 'basic');
    chi2    = xlsread(Options.weights_filepath, lc_ii, 'B2:Z2', 'basic');
    
    if all(chi2 == 0)
        Chi2q_ii        = chi2*ones(size(Chi1q_ii));
        Wq_ii           = interp1(chi1, W, chi1_ds);  
    else
        [Chi1, Chi2]    = meshgrid(chi1, chi2);
        Wq_ii           = interp2(Chi1, Chi2, W, Chi1q_ii, Chi2q_ii);
    end

    DS.Weight.(lc_ii).W     = lcw_ii*Wq_ii;
    DS.Weight.(lc_ii).Chi1  = Chi1q_ii;
    DS.Weight.(lc_ii).Chi2  = Chi2q_ii;
end 

% drop DS discretizations that are not relevant for the selected resistance
% model
DS = update_DS(DS, Options);

% =========================================================================
% PROBABILISTIC MODEL
% =========================================================================

% -------------------------------------------------------------------------
% RESISTANCE
% -------------------------------------------------------------------------

% theta_R factor in EC2 shear formula ~ resistance model uncertainty [-]
switch lower(resistance_model)
    case 'ec2_codified_2019'
        if consider_VRmin
            Prob.theta_R.mean   = 1.13688;
            Prob.theta_R.cov    = 0.23777;
%             Prob.theta_R.repr   = 1.0; % codified value
%             Prob.theta_R.repr   = 0.84679; % to obtain characteristics value (5%) for V_Rc with all relevant DSs
            Prob.theta_R.repr   = 0.84604; % to obtain characteristics value (5%) for V_Rc with the illustrative (reduced) set of DSs
        else
            Prob.theta_R.mean   = 1.13750;
            Prob.theta_R.cov    = 0.23760;
%             Prob.theta_R.repr   = 1.0; % codified value
            Prob.theta_R.repr   = 0.81776; % to obtain characteristics value (5%) for V_Rc
        end
        
        Prob.theta_R.std    = NaN;
        Prob.theta_R.dist   = 2;
        Prob.theta_R.P_repr = NaN;
        Prob.theta_R.gamma  = NaN;

    case 'ec2_pre_2021'
        % to be added
    case 'mc2010_level_ii_codified_2019'
        Prob.theta_R.mean   = 1.34439;
        Prob.theta_R.cov    = 0.19243;
%         Prob.theta_R.repr   = 0.96329; % to obtain characteristics value (5%) for V_Rc with all relevant DSs
%         Prob.theta_R.repr   = 1.07545; % to obtain characteristics value (5%) for V_Rc with all relevant DSs but only fc = [20, 40]
        Prob.theta_R.repr   = 1.07921; % to obtain characteristics value (5%) for V_Rc with the illustrative (reduced) set of DSs

        Prob.theta_R.std    = NaN;
        Prob.theta_R.dist   = 2;
        Prob.theta_R.P_repr = NaN;
        Prob.theta_R.gamma  = NaN;
    otherwise
        error(['Unknown resistance model:', resistance_model])
end

% The mean values are not relevant here, they will be overwritten by user
% specified values in main_calibr.m
% concrete compressive strength [MPa]
Prob.f_cc.mean      = 35;
Prob.f_cc.cov       = 0.15;
% old (2020)
% Prob.f_cc.cov       = 0.06;

Prob.f_cc.std       = NaN;
Prob.f_cc.dist      = 2;
Prob.f_cc.P_repr    = 0.05;
Prob.f_cc.gamma     = 1.5;
Prob.f_cc.repr2mean = @(x_repr, prob) repr2mean_fractile(x_repr, prob);
Prob.f_cc.mean2repr = @(x_mean, prob) mean2repr_fractile(x_mean, prob);

% effective depth [mm]
Prob.d.mean         = 200;
Prob.d.cov          = NaN;
Prob.d.std          = 10;
% old (2020)
% Prob.d.cov       = 5/200;

Prob.d.dist         = 1;
Prob.d.P_repr       = NaN;
Prob.d.gamma        = NaN;
Prob.d.repr2mean    = @(x_repr, prob) repr2mean_shift(x_repr, prob, 10);
Prob.d.mean2repr    = @(x_mean, prob) mean2repr_shift(x_mean, prob, 10);
% old (2020)
% Prob.d.repr2mean    = @(x_repr, prob) repr2mean_shift(x_repr, prob, 0);
% Prob.d.mean2repr    = @(x_mean, prob) mean2repr_shift(x_mean, prob, 0);

% width [mm]
Prob.b.mean         = 100;
Prob.b.cov          = NaN;
Prob.b.std          = 5;
% old (2020)
% Prob.b.cov       = 5/100;

Prob.b.dist         = 1;
Prob.b.P_repr       = NaN;
Prob.b.gamma        = NaN;
Prob.b.repr2mean    = @(x_repr, prob) repr2mean_shift(x_repr, prob, 0);
Prob.b.mean2repr    = @(x_mean, prob) mean2repr_shift(x_mean, prob, 0);

% % steel yield stress [MPa]
% Prob.f_sy.mean      = 440;
% Prob.f_sy.cov       = 0.06;
% Prob.f_sy.dist      = 2;
% Prob.f_sy.P_repr    = 0.045;
% Prob.f_sy.gamma     = 1.15;

% total rebar area [mm2]
Prob.Asl.mean       = 0.01*Prob.d.mean*Prob.b.mean;
Prob.Asl.cov        = 0.02;
Prob.Asl.dist       = 1;
Prob.Asl.P_repr     = NaN;
Prob.Asl.gamma      = NaN;
% this is more general as it works for asymmetric distributions as well
Prob.Asl.repr2mean  = @(x_repr, prob) repr2mean_shift(x_repr, prob, 0);
Prob.Asl.mean2repr  = @(x_mean, prob) mean2repr_shift(x_mean, prob, 0);

% size of the maximum sieve gird used for aggregates
% the mean (repr) value will be assigned during the optimization, comes
% from design scenario
Prob.d_lower.mean   = NaN;
Prob.d_lower.mean   = NaN;
Prob.d_lower.dist   = 0;

% a-d ratio
% the mean (repr) value will be assigned during the optimization, comes
% from design scenario
Prob.a_to_d_ratio.mean   = NaN;
Prob.a_to_d_ratio.repr   = NaN;
Prob.a_to_d_ratio.dist   = 0;

% -------------------------------------------------------------------------
% LOAD
% -------------------------------------------------------------------------

% Permanent load [kN]
Prob.G.mean         = 20;
Prob.G.cov          = 0.10;
Prob.G.std          = NaN;
Prob.G.dist         = 1;
Prob.G.P_repr       = 0.5;
Prob.G.gamma        = 1.35;
Prob.G.repr2mean    = @(x_repr, prob) repr2mean_fractile(x_repr, prob);
Prob.G.mean2repr    = @(x_mean, prob) mean2repr_fractile(x_mean, prob);

% Permanent load - effect model uncertainty [-]
Prob.K_G.mean       = 1.0;
Prob.K_G.cov        = 0.05;
Prob.K_G.std        = NaN;
Prob.K_G.dist       = 2;
Prob.K_G.P_repr     = NaN;
Prob.K_G.gamma      = NaN;

% Permanent load - ksi [-]
Prob.ksi.mean       = 0.85;
Prob.ksi.cov        = NaN;
Prob.ksi.dist       = 0;
Prob.ksi.P_repr     = NaN;
Prob.ksi.gamma      = NaN;

% Allocation of the variable loads 1 and 2, their model uncertainties and 
% psi-values (will be defined in inv_design.m)

% pseudo model [-]
pseudo.mean         = 0.0;
pseudo.cov          = NaN;
pseudo.dist         = 0;
pseudo.P_repr       = NaN;
pseudo.gamma        = 0.0;

Prob.Q1             = pseudo;
Prob.K_Q1           = pseudo;
Prob.psi01          = pseudo;
Prob.Q2             = pseudo;
Prob.K_Q2           = pseudo;
Prob.psi02          = pseudo;


% ACTION MODELS - Prob_actions

% Traffic load (annual maxima) - intensity [kN]
Prob_actions.T.mean        = 1.0;
Prob_actions.T.cov         = 0.075;
Prob_actions.T.P_repr      = 1-2.1257193880352965e-05; % to get 0.999 from the product of T*K_T;
% old (2020)
% Prob_actions.T.cov         = 0.0858;
% Prob_actions.T.P_repr      = 0.999;

Prob_actions.T.std         = NaN;
Prob_actions.T.dist        = 11;
Prob_actions.T.gamma       = 1.35;
Prob_actions.T.repr2mean   = @(x_repr, prob) repr2mean_fractile(x_repr, prob);
Prob_actions.T.mean2repr   = @(x_mean, prob) mean2repr_fractile(x_mean, prob);

% Traffic load - model uncertainty [-]
Prob_actions.K_T.mean      = 1.0;
Prob_actions.K_T.cov       = 0.142;
% old (2020)
% Prob_actions.K_T.cov       = 0.15;

Prob_actions.K_T.dist      = 1;
Prob_actions.K_T.P_repr    = NaN;
Prob_actions.K_T.gamma     = NaN;

% Traffic load - factor for combination value... EN1990 terminology [-]
Prob_actions.psi_T.mean     = 0.8;
Prob_actions.psi_T.cov      = NaN;
Prob_actions.psi_T.dist     = 0;
Prob_actions.psi_T.P_repr   = NaN;
Prob_actions.psi_T.gamma    = NaN;

% Snow load (annual maxima * ground2roof) - intensity
Prob_actions.S.mean        = NaN;
Prob_actions.S.cov         = NaN;
Prob_actions.S.dist        = 33;
Prob_actions.S.P_repr      = 0.97609;
% old (2020)
% Prob_actions.S.P_repr      = 0.97472;

Prob_actions.S.gamma       = 1.5;
Prob_actions.S.dist_ID     = 301;
Prob_actions.S.repr2mean   = @(x_repr, prob) repr2mean_fractile(x_repr, prob);
Prob_actions.S.mean2repr   = @(x_mean, prob) mean2repr_fractile(x_mean, prob);

% dummy test
% % Prob_actions.S.mean        = NaN;
% % Prob_actions.S.cov         = NaN;
% % Prob_actions.S.dist        = 33;
% % Prob_actions.S.P_repr      = 0.98;
% % Prob_actions.S.gamma       = 1.5;
% % Prob_actions.S.dist_ID     = 1001;

% % Prob_actions.S.mean        = 1;
% % Prob_actions.S.cov         = 0.4;
% % Prob_actions.S.dist        = 2;
% % Prob_actions.S.P_repr      = 0.98;
% % Prob_actions.S.gamma       = 1.5;

% Snow load - model uncertainty [-]
Prob_actions.K_S.mean      = 1;
Prob_actions.K_S.cov       = 0.10;
Prob_actions.K_S.dist      = 2;
Prob_actions.K_S.P_repr    = NaN;
Prob_actions.K_S.gamma     = NaN;

% Snow load - factor for combination value... EN1990 terminology [-]
Prob_actions.psi_S.mean     = 0.5;
Prob_actions.psi_S.cov      = NaN;
Prob_actions.psi_S.dist     = 0;
Prob_actions.psi_S.P_repr   = NaN;
Prob_actions.psi_S.gamma    = NaN;

% Wind load (annual maxima*all the rest but model unc) - intensity
Prob_actions.W.mean        = NaN;
Prob_actions.W.cov         = NaN;
Prob_actions.W.dist        = 33;
Prob_actions.W.P_repr      = 0.99037;
Prob_actions.W.gamma       = 1.50;
Prob_actions.W.dist_ID     = 201;
Prob_actions.W.repr2mean   = @(x_repr, prob) repr2mean_fractile(x_repr, prob);
Prob_actions.W.mean2repr   = @(x_mean, prob) mean2repr_fractile(x_mean, prob);

% Wind load - time-invariant component with model uncertainty [-]
Prob_actions.K_W.mean      = 1.0;
Prob_actions.K_W.cov       = 0.10;
Prob_actions.K_W.dist      = 2;
Prob_actions.K_W.P_repr    = NaN;
Prob_actions.K_W.gamma     = NaN;

% Wind load - factor for combination value... EN1990 terminology [-]
Prob_actions.psi_W.mean     = 0.6;
Prob_actions.psi_W.cov      = NaN;
Prob_actions.psi_W.dist     = 0;
Prob_actions.psi_W.P_repr   = NaN;
Prob_actions.psi_W.gamma    = NaN;

% Imposed load (annual maxima) - intensity [kN]
Prob_actions.I.mean        = 1.0;
Prob_actions.I.cov         = 0.53;
Prob_actions.I.dist        = 11;
Prob_actions.I.P_repr      = 0.98;
Prob_actions.I.gamma       = 1.50;
Prob_actions.I.repr2mean   = @(x_repr, prob) repr2mean_fractile(x_repr, prob);
Prob_actions.I.mean2repr   = @(x_mean, prob) mean2repr_fractile(x_mean, prob);

% Imposed load - model uncertainty [-]
Prob_actions.K_I.mean      = 1.0;
Prob_actions.K_I.cov       = 0.10;
Prob_actions.K_I.dist      = 2;
Prob_actions.K_I.P_repr    = NaN;
Prob_actions.K_I.gamma     = NaN;

% Imposed load - factor for combination value... EN1990 terminology [-]
Prob_actions.psi_I.mean     = 0.7;
Prob_actions.psi_I.cov      = NaN;
Prob_actions.psi_I.dist     = 0;
Prob_actions.psi_I.P_repr   = NaN;
Prob_actions.psi_I.gamma    = NaN;

% Load effect model uncertainty [-]
Prob_actions.K_E.mean      = 1.0;
Prob_actions.K_E.cov       = 0.10;
% old (2020)
% Prob_actions.K_E.cov       = 0.00010;
Prob_actions.K_E.dist      = 2;
Prob_actions.K_E.P_repr    = NaN;
Prob_actions.K_E.gamma     = NaN;

Prob.K_E = Prob_actions.K_E;

end
