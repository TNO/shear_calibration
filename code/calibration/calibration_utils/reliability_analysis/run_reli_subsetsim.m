% Reliability analysis FERUM: subset simulation
%
%SYNOPSYS
% [beta, simulationresults, probdata] = RUN_RELI_subsetsim(Prob, Options)
%
%

%NOTES:

function [beta, subsetsimulationresults, probdata] = run_reli_subsetsim(Prob, Options)

% =========================================================================
% INPUT DATA
% =========================================================================

resistance_model            = Options.resistance_model;
load_combination            = Options.load_combination;

rv_order                    = {'C', 'f_cc', 'd', 'b', 'Asl', 'G', 'K_G', 'ksi', 'Q1', 'K_Q1', 'psi01', 'Q2', 'K_Q2', 'psi02'};    
Prob                        = orderfields(Prob, rv_order);

var_names                   = fieldnames(Prob);
n_var                       = length(var_names);

% =========================================================================
% FERUM
% =========================================================================

% -------------------------------------------------------------------------
%  DATA FIELDS IN 'PROBDATA'
% -------------------------------------------------------------------------

% Names of random variables. Default names are 'x1', 'x2', ..., if not explicitely defined.

probdata.name               = var_names;

% Marginal distributions for each random variable
% probdata.marg =  [ (type) (mean) (stdv) (startpoint) (p1) (p2) (p3) (p4) (input_type); ... ];

% s_g's distribution is defined with distribution parameters, input_type = 1!
marg                        = nan(n_var, 9);
for ii = 1:n_var
    X          = Prob.(var_names{ii});
    marg(ii,:) = [X.dist,  X.mean,  X.mean*X.cov,  X.mean,  NaN,  NaN,  NaN,  NaN, 0];
end        

% .........................................................................
% Auxiliary variables
% .........................................................................
% resistance model ID
resi_model                      = translate_model(resistance_model);
marg(n_var+1,:)                 = [0,  resi_model,  0,  resi_model,  NaN,  NaN,  NaN,  NaN, 0];

% resistance model ID
load_comb                       = translate_model(load_combination);
marg(n_var+2,:)                 = [0,  load_comb,  0,  load_comb,  NaN,  NaN,  NaN,  NaN, 0];

probdata.name                   = {probdata.name{:}, 'resi_model', 'load_comb'}'; 
n_var                           = n_var + 2;

probdata.marg                   = marg;

% Correlation matrix
probdata.correlation            = eye(n_var);
%rv_order: {'C', 'f_cc', 'd', 'b', 'Asl', 'G', 'K_G', 'ksi', 'Q1', 'K_Q1', 'psi01', 'Q2', 'K_Q2', 'psi02'};    
%probdata.correlation(3,6) = 0.99; probdata.correlation(6,3) = probdata.correlation(3,6);
%probdata.correlation(4,6) = 0.70; probdata.correlation(6,4) = probdata.correlation(4,6);

probdata.transf_type            = 3;
probdata.Ro_method              = 1;
probdata.flag_sens              = 0;

% -------------------------------------------------------------------------
%  DATA FIELDS IN 'ANALYSISOPT'
% -------------------------------------------------------------------------

analysisopt.echo_flag            = 0;

analysisopt.multi_proc           = 1;        % 1: block_size g-calls sent simultaneously
%    - gfunbasic.m is used and a vectorized version of gfundata.expression is available.
%      The number of g-calls sent simultaneously (block_size) depends on the memory
%      available on the computer running FERUM.
%    - gfunxxx.m user-specific g-function is used and able to handle block_size computations
%      sent simultaneously, on a cluster of PCs or any other multiprocessor computer platform.
% 0: g-calls sent sequentially
analysisopt.block_size           = 1e4;     % Number of g-calls to be sent simultaneously

% Simulation analysis (MC,IS,DS,SS) and distribution analysis options
analysisopt.num_sim              = 5e4   ;   % Number of samples (MC,IS), number of samples per subset step (SS) or number of directions (DS)
analysisopt.rand_generator       = 0;        % 0: default rand matlab function, 1: Mersenne Twister (to be preferred)
%disp('random generator Mersenne Twister is not working')

% Simulation analysis (MC, IS) and distribution analysis options
%analysisopt.sim_point            = 'dspt';  % 'dspt': design point, 'origin': origin in standard normal space (simulation analysis)
analysisopt.stdv_sim             = 1;        % Standard deviation of sampling distribution in simulation analysis

% Simulation analysis (MC, IS)
analysisopt.target_cov           = 0.01;   % Target coefficient of variation for failure probability
%analysisopt.lowRAM               = 0;        % 1: memory savings allowed, 0: no memory savings allowed

% Subset Simulation (SS) analysis options
analysisopt.width                = 2;        % Width of the proposal uniform pdfs
analysisopt.pf_target            = 0.1;      % Target probability for each subset step
analysisopt.flag_cov_pf_bounds   = 1;        % 1: calculate upper and lower bounds of the coefficient of variation of pf
                                             % 0: no calculation
analysisopt.ss_restart_from_step = -inf;     % i>=0 : restart from step i
                                             % -inf : all steps, no record (default)
                                             % -1 : all steps, record all
analysisopt.flag_plot            = 0;        % 1: plots at each step (2 r.v. examples only)
                                             % 0: no plots
analysisopt.flag_plot_gen        = 0;        % 1: intermediate plots for each MCMC chain (2 r.v. examples only)
                                             % 0: no plots

% -------------------------------------------------------------------------
%  DATA FIELDS IN 'GFUNDATA' (one structure per gfun) 
% -------------------------------------------------------------------------

% Type of limit-state function evaluator:
% 'basic': the limit-state function is defined by means of an analytical expression or a Matlab m-function,
%          using gfundata(lsf).expression. The function gfun.m calls gfunbasic.m, which evaluates gfundata(lsf).expression.
% 'xxx':   the limit-state function evaluation requires a call to an external code.  The function gfun.m calls gfunxxx.m,
%          which evaluates gfundata(lsf).expression where gext variable is a result of the external code.
gfundata(1).evaluator           = 'basic';
gfundata(1).type                = 'expression';   % Do not change this field!

% Expression of the limit-state function:


% gfundata(1).expression = 'gfun_diana(fcc, Ec)';
arg                             = sprintf('%s,', probdata.name{:});
arg(end)                        = [];

gfundata(1).expression          = ['g_fun(',arg, ')'];


gfundata(1).thetag              =  [];
gfundata(1).thetagname          = { };


% Flag for computation of sensitivities w.r.t. thetag parameters of the limit-state function
% 1: all sensitivities assessed, 0: no sensitivities assessment
gfundata(1).flag_sens           = 0;

% -------------------------------------------------------------------------
%  DATA FIELDS IN 'FEMODEL'
% -------------------------------------------------------------------------

femodel = [];

% -------------------------------------------------------------------------
%  DATA FIELDS IN 'RANDOMFIELD'
% -------------------------------------------------------------------------

randomfield = [];

% -------------------------------------------------------------------------
% Subset simulation ANALYSIS
% -------------------------------------------------------------------------

% This function updates probdata and gfundata before any analysis (must be run only once)
[probdata, gfundata, analysisopt] = update_data(1, probdata, analysisopt, gfundata, []);

% This function completely determines and updates parameters, mean and standard deviation associated with the distribution of each random variable
probdata.marg = distribution_parameter(probdata.marg);

% Subset simulation:
%tic;
[ subsetsimulationresults, probdata ] = subset_simulation(1,probdata,analysisopt,gfundata,femodel,randomfield);
%toc

beta = subsetsimulationresults.beta;

