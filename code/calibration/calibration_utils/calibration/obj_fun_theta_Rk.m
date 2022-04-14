%
%

function [fval, Results, DS] = obj_fun_theta_Rk(x, Prob, Prob_actions, DS, Options)

% -------------------------------------------------------------------------
% Initialize
% -------------------------------------------------------------------------
P_repr_target           = Options.P_repr_target;
Prob.theta_R.repr       = x;

keep_fields = {'theta_R', 'f_cc', 'd', 'b', 'Asl', 'd_lower', 'a_to_d_ratio'};
all_fields = fieldnames(Prob);
remove_fields = setdiff(all_fields, keep_fields);

% -------------------------------------------------------------------------
% Semi-probabilistic design
% -------------------------------------------------------------------------
% semi-probabilistic design
free_par                = 1.0;
[Prob, DS]              = gen_DS(free_par, Prob, Prob_actions, DS, Options);

n_ds                    = size(DS.p_ds_all, 1);

% -------------------------------------------------------------------------
% Reliability analysis
% -------------------------------------------------------------------------
beta                    = nan(n_ds,1);
alphas                  = cell(n_ds,1);
converged               = false(n_ds,1);
parfor ii = 1:n_ds
% for ii = 1:n_ds
    Prob_ii                = Prob(ii);
    % keep only the relevant RVs
    Prob_ii = rmfield(Prob_ii, remove_fields);

    [beta(ii),formresults] = run_reli_theta_Rk(Prob_ii, Options);
    alphas{ii}             = formresults.alpha;
    converged(ii)          = formresults.converged;

    if formresults.converged == 0
        error('FORM did not converge.')
    end
end

% -------------------------------------------------------------------------
% Objective function
% -------------------------------------------------------------------------
% The commented out formulation is correct until the objective goes down
% to zero, if that is no longer the case then `weighted_integral` 
% should be used.
% P_repr = normcdf(-beta);
% weights = DS.weights_combis;
% fval  = trapz(weights.*(P_repr_target - P_repr).^2);
Options.beta_target = -norminv(P_repr_target);
Options.objective_function = "squared";
fval = weighted_integral(beta, DS, Options);

% -------------------------------------------------------------------------
% Collect results
% -------------------------------------------------------------------------
Results.alphas = alphas;
Results.beta            = beta;
Results.converged       = converged;

end