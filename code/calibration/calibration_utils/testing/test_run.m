clearvars
close all
clc

Options.resistance_model    = 'ec2';
Options.load_combination    = 'ec2-simple';
%Options.load_combination    = 'ec2-advanced';
Options.beta_target         = 3.8;
Options.verbose             = 2;

%--------------------------------------------------------------------------
% INPUT DATA
%--------------------------------------------------------------------------
free_par                    = 1.2904;
[Prob, DS]                  = get_input(Options);
Prob                        = update_Prob(Prob);

make_table(Prob)

%--------------------------------------------------------------------------
% TEST THE OBJ FUN
%--------------------------------------------------------------------------
% [fval, Results]             = obj_fun(free_par, Prob, DS, Options)

%--------------------------------------------------------------------------
% SEMI-PROB DESIGN
%--------------------------------------------------------------------------

[Prob, DS]                  = gen_DS(free_par, Prob, DS, Options);

n_ds                        = length(Prob);
%--------------------------------------------------------------------------
% RELIABILITY ANALYSIS
%--------------------------------------------------------------------------
dsnumber = 1
[beta1, formresults1, probdata1]          = run_reli(Prob(dsnumber), Options);
%[beta, formresults, probdata]             = run_reli_cobyla(Prob(dsnumber), Options);
%[beta, simulationresults, probdata]       = run_reli_cMC(Prob(dsnumber), Options);
[beta, dirsimulationresults, probdata]    = run_reli_dirsim(Prob(dsnumber), Options);
[beta, subsetsimulationresults, probdata] = run_reli_subsetsim(Prob(dsnumber), Options);

formresults1
% formresults
% simulationresults
dirsimulationresults
subsetsimulationresults

% Beta = nan(n_ds,1);
% for ii = 1:n_ds
%     [Beta(ii), formresults, probdata] = run_reli(Prob(ii), Options);
% end
% 
% formresults
% 
% plot(Beta, '-o')
