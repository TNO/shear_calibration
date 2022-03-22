%
%

function [fval, Results, DS] = obj_fun_Ck(x, Prob, Prob_actions, DS, Options)

% -------------------------------------------------------------------------
% Initialize
% -------------------------------------------------------------------------
verbose                 = Options.verbose;
P_repr_target           = Options.P_repr_target;
Prob.C.in_standardized_equation = x;

keep_fields = {'C', 'f_cc', 'd', 'b', 'Asl'};
all_fields = fieldnames(Prob);
remove_fields = setdiff(all_fields, keep_fields);
% -------------------------------------------------------------------------
% Semi-probabilistic design
% -------------------------------------------------------------------------
% semi-probabilistic design
free_par                = 1.0;
[Prob, DS]              = gen_DS(free_par, Prob, Prob_actions, DS, Options);

n_ds                    = size(DS.combis, 1);

% -------------------------------------------------------------------------
% Reliability analysis
% -------------------------------------------------------------------------
beta                    = nan(n_ds,1);
alphas                  = cell(n_ds,1);
flag_form               = false(n_ds,1);
parfor ii = 1:n_ds
% for ii = 1:n_ds
    Prob_ii                = Prob(ii);
    % keep only the relevant RVs
    Prob_ii = rmfield(Prob_ii, remove_fields);

    [beta(ii),formresults] = run_reli_Ck(Prob_ii, Options);
    alphas{ii}             = formresults.alpha;
    flag_form(ii)          = formresults.flag;
end

% -------------------------------------------------------------------------
% Objective function
% -------------------------------------------------------------------------
P_repr = normcdf(-beta);
weights = DS.weights_combis;
fval  = trapz(weights.*(P_repr_target - P_repr).^2);

% -------------------------------------------------------------------------
% Collect results
% -------------------------------------------------------------------------
if verbose == 0
   Results.alphas = alphas;
end
Results.beta            = beta;
Results.flag_form       = flag_form;

end