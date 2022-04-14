% Reliability analysis FERUM::FORM
%
%SYNOPSYS
% [beta, formresults, probdata] = RUN_RELI(Prob, Options)
%
%

%NOTES:
% - vectorized g_fun calls

function [beta, reli_results, probdata] = run_reli(Prob, Options, method, formresults)

if nargin < 4 && strcmp(method, 'is')
    error('For Importance sampling (is) `formresults` should be provided.')
end

% =========================================================================
% INPUT DATA
% =========================================================================

resistance_model            = Options.resistance_model;
load_combination            = Options.load_combination;
consider_VRmin              = Options.consider_VRmin;

rv_order                    = {
    'theta_R', 'f_cc', 'd', 'b', 'Asl', 'd_lower', 'a_to_d_ratio', ...
    'G', 'K_G', 'ksi', 'Q1', 'K_Q1', 'psi01', 'Q2', 'K_Q2', 'psi02', 'K_E'
};    
Prob                        = orderfields(Prob, rv_order);

var_names                   = fieldnames(Prob);
n_var                       = length(var_names);

% e.g. if VRmin and VRbase are both considered in the EC2 model than that
% more than one branch; the reliability method is adjusted based on this
% boolean variable (see below the limit state functiond definition).
more_than_one_resistance_branch = ...
    any(strcmpi(resistance_model, {'ec2_codified_2019'})) && consider_VRmin == true;

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

    x_dist = X.dist;
    x_mean = X.mean;
    x_cov = X.cov;
    x_std = X.std;

    switch x_dist
        case {1, 2, 3, 11}
            if (~isnan(x_cov) && ~isnan(x_std)) || (isnan(x_cov) && isnan(x_std))
                error('Only one of the `X.cov` `X.std` pairs must differ from NaN.')
            end
        
        if isnan(x_std)
            x_std = x_cov * x_mean;
            X.std = x_std;
        end

        case {0, 33}
            % do nothing
        otherwise
            error(['Unkown or not implemented distributiont type (dist): ', num2str(x_dist)])
    end

    marg(ii,:) = [X.dist,  X.mean,  X.std,  X.mean,  X.dist_ID,  X.shift,  X.scale,  NaN, 0];
end        

% .........................................................................
% Auxiliary variables
% .........................................................................
% resistance model ID
resi_model                      = translate_model(resistance_model);
marg(n_var+1,:)                 = [0,  resi_model,  0,  resi_model,  NaN,  NaN,  NaN,  NaN, 0];

% load combination ID
load_comb                       = translate_model(load_combination);
marg(n_var+2,:)                 = [0,  load_comb,  0,  load_comb,  NaN,  NaN,  NaN,  NaN, 0];

% option for the resistance formula
% marg(n_var+3,:)                 = [0,  consider_VRmin,  0,  consider_VRmin,  NaN,  NaN,  NaN,  NaN, 0];

probdata.name                   = [probdata.name(:)', {'resi_model'}, {'load_comb'}]'; 
n_var                           = n_var + 3;

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

% FORM analysis options
analysisopt.i_max                = 500;      % Maximum number of iterations allowed in the search algorithm
analysisopt.e1                   = 0.01;    % Tolerance on how close design point is to limit-state surface
analysisopt.e2                   = 0.01;    % Tolerance on how accurately the gradient points towards the origin
analysisopt.step_code            = 0;        % 0: step size by Armijo rule, otherwise: given value is the step size
analysisopt.Recorded_u           = 1;        % 0: u-vector not recorded at all iterations, 1: u-vector recorded at all iterations
analysisopt.Recorded_x           = 1;        % 0: x-vector not recorded at all iterations, 1: x-vector recorded at all iterations

% FORM, SORM analysis options
analysisopt.grad_flag            = 'ffd';    % 'ddm': direct differentiation, 'ffd': forward finite difference, 'ffd_ddm' mixed forward finite difference and direct differentiation, NaN in dgdq indicates that ffd is needed for the estimation!
% analysisopt.grad_flag            = 'ffd_ddm';    % 'ddm': direct differentiation, 'ffd': forward finite difference, 'ffd_ddm' mixed forward finite difference and direct differentiation, NaN in dgdq indicates that ffd is needed for the estimation!
analysisopt.ffdpara              = 1000;     % Parameter for computation of FFD estimates of gradients - Perturbation = stdv/analysisopt.ffdpara;
% Recommended values: 1000 for basic limit-state functions, 50 for FE-based limit-state functions
analysisopt.ffdpara_thetag       = 1000;     % Parameter for computation of FFD estimates of dbeta_dthetag
% perturbation = thetag/analysisopt.ffdpara_thetag if thetag ~= 0 or 1/analysisopt.ffdpara_thetag if thetag == 0;
% Recommended values: 1000 for basic limit-state functions, 100 for FE-based limit-state functions

%analysisopt.form_solver             = 'cobyla';

% Simulation analysis (MC,IS,DS,SS) and distribution analysis options
analysisopt.num_sim              = 1e5;      % Number of samples (MC,IS), number of samples per subset step (SS) or number of directions (DS)
analysisopt.rand_generator       = 1;        % 0: default rand matlab function, 1: Mersenne Twister (to be preferred)

% Simulation analysis (MC, IS) and distribution analysis options
analysisopt.sim_point            = 'dspt';   % 'dspt': design point, 'origin': origin in standard normal space (simulation analysis)
analysisopt.stdv_sim             = 1;        % Standard deviation of sampling distribution in simulation analysis

% Subset Simulation (SS) analysis options
analysisopt.width                = 2;        % Width of the proposal uniform pdfs
analysisopt.pf_target            = 0.1;      % Target probability for each subset step
analysisopt.flag_cov_pf_bounds   = 1;        % 1: calculate upper and lower bounds of the coefficient of variation of pf
                                             % 0: no calculation
analysisopt.ss_restart_from_step = -inf;     % i>=0 : restart from step i
                                             % -inf : all steps, no record (default)
                                             % -1 : all steps, record all

% Simulation analysis (MC, IS)
analysisopt.target_cov           = 0.01;   % Target coefficient of variation for failure probability
analysisopt.lowRAM               = 0;        % 1: memory savings allowed, 0: no memory savings allowed

% -------------------------------------------------------------------------
%  DATA FIELDS IN 'GFUNDATA' (one structure per gfun) 
% -------------------------------------------------------------------------

% Type of limit-state function evaluator:
% 'basic': the limit-state function is defined by means of an analytical expression or a Matlab m-function,
%          using gfundata(lsf).expression. The function gfun.m calls gfunbasic.m, which evaluates gfundata(lsf).expression.
% 'xxx':   the limit-state function evaluation requires a call to an external code.  The function gfun.m calls gfunxxx.m,
%          which evaluates gfundata(lsf).expression where gext variable is a result of the external code.

if more_than_one_resistance_branch == true
    % .....................................................................
    % LSF1: only VRbase is considered
    % .....................................................................
    consider_VRbase_lsf = true;
    consider_VRmin_lsf = false;
    
    gfundata(1).evaluator           = 'basic';
    gfundata(1).type                = 'expression';   % Do not change this field!
    
    arg                             = sprintf('%s,', probdata.name{:});
    
    gfundata(1).expression          = [
        'g_fun(', arg, num2str(consider_VRmin_lsf), ',', ...
        num2str(consider_VRbase_lsf) ,')'];
    
    gfundata(1).thetag              = [];
    gfundata(1).thetagname          = {};
    % Flag for computation of sensitivities w.r.t. thetag parameters of 
    % the limit-state function
    % 1: all sensitivities assessed, 0: no sensitivities assessment
    gfundata(1).flag_sens           = 0;
    
    % .....................................................................
    % LSF2: only VRmin is not considered
    % .....................................................................
    consider_VRbase_lsf = false;
    consider_VRmin_lsf = true;
    
    gfundata(2) = gfundata(1);
    gfundata(2).expression          = [
        'g_fun(', arg, num2str(consider_VRmin_lsf), ',', ...
        num2str(consider_VRbase_lsf) ,')'];

else
    consider_VRbase_lsf             = true;

    gfundata(1).evaluator           = 'basic';
    gfundata(1).type                = 'expression';   % Do not change this field!
    
    arg                             = sprintf('%s,', probdata.name{:});
    
    gfundata(1).expression          = [
        'g_fun(', arg, num2str(consider_VRmin), ',', ...
        num2str(consider_VRbase_lsf) ,')'];

    gfundata(1).thetag              = [];
    gfundata(1).thetagname          = {};
    % Flag for computation of sensitivities w.r.t. thetag parameters of 
    % the limit-state function
    % 1: all sensitivities assessed, 0: no sensitivities assessment
    gfundata(1).flag_sens           = 0;
end

% -------------------------------------------------------------------------
%  DATA FIELDS IN 'FEMODEL'
% -------------------------------------------------------------------------

femodel = [];

% -------------------------------------------------------------------------
%  DATA FIELDS IN 'RANDOMFIELD'
% -------------------------------------------------------------------------

randomfield = [];

% -------------------------------------------------------------------------
% RELIABILITY ANALYSIS
% -------------------------------------------------------------------------

switch method
    case 'form'
        % FORM analysis
        if more_than_one_resistance_branch == true
            % .............................................................
            % LSF1
            % .............................................................
            % This function updates probdata and gfundata before any analysis (must be run only once)
            [probdata_1, gfundata_1, analysisopt_1] = update_data(1, probdata, analysisopt, gfundata(1), []);

            % This function completely determines and updates parameters, mean and standard deviation associated with the distribution of each random variable
            probdata_1.marg = distribution_parameter(probdata_1.marg);

            t_start = tic;
            [formresults_1, probdata_1] = form( ...
                1, probdata_1, analysisopt_1, gfundata_1, femodel, randomfield);
            t_1_elapsed_sec = toc(t_start);
            
            % .............................................................
            % LSF2
            % .............................................................
            % This function updates probdata and gfundata before any analysis (must be run only once)
            [probdata_2, gfundata_2, analysisopt_2] = update_data(1, probdata, analysisopt, gfundata(2), []);

            % This function completely determines and updates parameters, mean and standard deviation associated with the distribution of each random variable
            probdata_2.marg = distribution_parameter(probdata_2.marg);

            t_start = tic;
            [formresults_2, probdata_2] = form(1, probdata_2, analysisopt_2, gfundata_2, femodel, randomfield);
            t_2_elapsed_sec = toc(t_start);

            t_elapsed_sec = t_1_elapsed_sec + t_2_elapsed_sec;

            % .............................................................
            % Parallel system
            % .............................................................
            beta_1 = formresults_1.beta;
            beta_2 = formresults_2.beta;
            alpha_1 = formresults_1.alpha;
            alpha_2 = formresults_2.alpha;
            
            % correlation coefficient between the two failures (limit states)
            rho = alpha_1' * alpha_2;

            Pf_s = mvncdf(-[beta_1, beta_2], [0, 0], [1, rho; rho, 1]);
            beta_s = norminv(1 - Pf_s);

            % return the `formresults` of the limit state with the larger
            % (governing) reliability index: 
            %  - keep its alpha value as an approximation;
            %  - replace its beta (Pf) with the system beta (Pf)
            %  - replace its convergence flag with the combined convergence
            %    of both FORM analyses
            if beta_1 > beta_2
                formresults = formresults_1;
                probdata = probdata_1;
            else
                formresults = formresults_2;
                probdata = probdata_2;
            end

            formresults.beta = beta_s;
            formresults.pf1 = Pf_s;

            % convergence flag
            if formresults_1.iter < analysisopt.i_max && formresults_2.iter < analysisopt.i_max
                formresults.converged = true;
            else
                formresults.converged = false;
            end

        else
            % This function updates probdata and gfundata before any analysis (must be run only once)
            [probdata, gfundata, analysisopt] = update_data(1, probdata, analysisopt, gfundata, []);

            % This function completely determines and updates parameters, mean and standard deviation associated with the distribution of each random variable
            probdata.marg = distribution_parameter(probdata.marg);

            t_start = tic;
            [formresults, probdata] = form(1, probdata, analysisopt, gfundata, femodel, randomfield);
            t_elapsed_sec = toc(t_start);
            
            % convergence flag
            if formresults.iter < analysisopt.i_max
                formresults.converged = true;
            else
                formresults.converged = false;
            end
        end
        
        % collect results
        formresults.t_elapsed = t_elapsed_sec;
        reli_results = formresults;
        beta = formresults.beta;

    case 'is'
        % .................................................................
        % Importance sampling with FORM dspt
        % .................................................................

        % This function updates probdata and gfundata before any analysis (must be run only once)
        [probdata, gfundata, analysisopt] = update_data(1, probdata, analysisopt, gfundata, []);

        % This function completely determines and updates parameters, mean and standard deviation associated with the distribution of each random variable
        probdata.marg = distribution_parameter(probdata.marg);

        analysisopt.formresults = formresults;
        t_start = tic;
        [simulationresults, probdata] = simulation_single_dspt(1, probdata, analysisopt, gfundata, femodel, randomfield);
        t_elapsed_sec = toc(t_start);
        
        % convergence flag
        if simulationresults.cov_pf < analysisopt.target_cov
            simulationresults.converged = true;
        else
            simulationresults.converged = false;
        end

        simulationresults.t_elapsed = t_elapsed_sec;
        reli_results = simulationresults;
        beta = simulationresults.beta;

    case 'subset'
        % .................................................................
        % Subset simulation
        % .................................................................

        % This function updates probdata and gfundata before any analysis (must be run only once)
        [probdata, gfundata, analysisopt] = update_data(1, probdata, analysisopt, gfundata, []);

        % This function completely determines and updates parameters, mean and standard deviation associated with the distribution of each random variable
        probdata.marg = distribution_parameter(probdata.marg);

        t_start = tic;
        [subsetsimulationresults, probdata] = subset_simulation(1, probdata, analysisopt, gfundata, femodel, randomfield);
        t_elapsed_sec = toc(t_start);
        
        subsetsimulationresults.t_elapsed = t_elapsed_sec;
        reli_results = subsetsimulationresults;
        beta = subsetsimulationresults.beta;
    otherwise
        error(['Unknown reliability method: ', method])
end


