% ?j?m
%
%
%NOTES:
% * UNITS should be carefully checked because of the messed up resistance
%   formulas!!
% * assuming that ksi, psi01, and psi02 are related to the simultaneous occurance of actions 
%
%
%TODO
% * DS-Prob consistency check! and warning!

clearvars
close all
clc

to_path = {'calibration_utils\', '..\shared_utils\resistance_models\'};
cellfun(@(x) addpath(genpath(x)), to_path)

%--------------------------------------------------------------------------
% OPTIONS
%--------------------------------------------------------------------------
% Shear resistance model/formula
% 'ec2_codified_2019', 'ec2_new', 'ec2_proposed_tg4_2016', 'ec2_proposed_yuguang_2019'
% 'mc2010_level_ii_codified_2019', 'mc2010_new'
Options.resistance_model    = 'ec2_codified_2019';

Options.consider_VRmin      = true;
% Options.consider_VRmin      = false;

% .........................................................................
% Only to get the design scenarios; not used in the reli calculations
% Load combination rule/formula: 
% 'ec0_simple', 'ec0_advanced'
Options.load_combination    = 'ec0_simple';
% Options.load_combination    = 'ec0_advanced';

% Variable load sets (#load comb would be more descriptive)
% 'traffic', 'snow-wind', 'snow-imposed', 'wind-imposed'
Options.load_combs          = {'traffic', 'snow_wind', 'snow_imposed', 'wind_imposed'};
% Options.load_combs          = {'snow_wind'}; 
% Options.load_combs          = {'traffic', 'wind_imposed'};
% Options.load_combs          = {'traffic'};

Options.load_comb_weights   = [1, 1, 1, 1];

% relative to the location of this file
data_dir                    = '../../data/';
% Weights for design scenarios (comment or uncomment)
Options.weights_filepath    = fullfile(data_dir, 'load_comb_prevalence_weights.xlsx');

Options.K_FI_repr           = 1.0;
% .........................................................................

% Target reliability
Options.P_repr_target       = 0.05;

Options.verbose             = 1;

%--------------------------------------------------------------------------
% GET INPUT DATA
%--------------------------------------------------------------------------
[Prob, Prob_actions, DS]    = get_input(Options);
Prob                        = update_Prob(Prob, Options.verbose);
Prob_actions                = update_Prob(Prob_actions, Options.verbose);

%--------------------------------------------------------------------------
% CALIBRATE
%--------------------------------------------------------------------------
tic
[calibr_par, objfun_val, exitflag] = calibrate_theta_Rk(Prob, Prob_actions, DS, Options);
toc
disp(['calibrated theta_R_repr (P=',...
    sprintf('%.2f', Options.P_repr_target), '): ',...
    sprintf('%.5f', calibr_par)])

%--------------------------------------------------------------------------
% CLEAN UP
%--------------------------------------------------------------------------
% Remove from the path: for safety as the functions added to the path are
% available from everywhere
cellfun(@(x) rmpath(genpath(x)), to_path)