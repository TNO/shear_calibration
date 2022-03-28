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
% Options.resistance_model    = 'ec2_codified_2019';
% Options.resistance_model    = 'ec2_pre_2021';
Options.resistance_model    = 'mc2010_level_ii_codified_2019';

Options.consider_VRmin      = true;
% Options.consider_VRmin      = false;

% Load combination rule/formula
% 'ec0_simple', 'ec0_advanced'
Options.load_combination    = 'ec0_simple';
% Options.load_combination    = 'ec0_advanced';

% Variable load sets (#load comb would be more descriptive
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

% Target reliability
% Options.beta_target         = 4.2;
Options.beta_target         = 4.7;
% Options.beta_target         = 5.2;

% EN 1990, Table B3; action scaler in semi-probabilistic design
% RC1: beta_target = 4.2; K_FI_repr = 0.9
% RC2: beta_target = 4.7; K_FI_repr = 1.0
% RC3: beta_target = 5.2; K_FI_repr = 1.1
% Options.K_FI_repr           = 0.9;
Options.K_FI_repr           = 1.0;
% Options.K_FI_repr           = 1.1;

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
[calibr_par, objfun_val, exitflag] = calibrate(Prob, Prob_actions, DS, Options);
toc
disp('calibr_par')
disp(calibr_par)

Options.verbose             = 0;
[~, Results, DS]            = obj_fun(calibr_par, Prob, Prob_actions, DS, Options);

Results.calibr_par          = calibr_par;

%--------------------------------------------------------------------------
% SAVE
%--------------------------------------------------------------------------
ID                          = datestr(now,'YYYY-mmm-dd_HH.MM.SS');

save(['results\calibration_run_', ID, '.mat'], 'Results', 'Options', 'DS', 'Prob', 'Prob_actions')

% save the results for visualization
output = {'beta', 'alpha'};

for i = 1:length(output)
    fname  = ['calibration_run_', ID];
    save_results_for_visu(fname, output{i})
end

% quick diagnostic plot
quick_visual_results(Results, Options.beta_target, sprintf('calibrated, gamma_C=%.3f', calibr_par))

% .........................................................................
% reference with current partial factor
calibr_par                  = 1.5;
ID                          = [ID, '_base_', sprintf('%.1f', calibr_par)];
Options.verbose             = 0;
[~, Results, DS]            = obj_fun(calibr_par, Prob, Prob_actions, DS, Options);

Results.calibr_par          = calibr_par;
save(['results\calibration_run_', ID, '.mat'], 'Results', 'Options', 'DS', 'Prob', 'Prob_actions')

% save the results for visualization
output = {'beta', 'alpha'};

for i = 1:length(output)
    fname  = ['calibration_run_', ID];
    save_results_for_visu(fname, output{i})
end

% quick diagnostic plot
quick_visual_results(Results, Options.beta_target, sprintf('reference, gamma_C=%.3f', calibr_par))

%--------------------------------------------------------------------------
% CLEAN UP
%--------------------------------------------------------------------------
% Remove from the path: for safety as the functions added to the path are
% available from everywhere
cellfun(@(x) rmpath(genpath(x)), to_path)