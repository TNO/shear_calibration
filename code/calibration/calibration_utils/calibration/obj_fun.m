%TODO
% * collect all formresults into an array of structure
%
%

function [fval, Results, DS] = obj_fun(x, Prob, Prob_actions, DS, Options)

% -------------------------------------------------------------------------
% Initialize
% -------------------------------------------------------------------------
verbose                 = Options.verbose;
beta_t                  = Options.beta_target;

% -------------------------------------------------------------------------
% Semi-probabilistic design
% -------------------------------------------------------------------------
% semi-probabilistic design
free_par                = x;
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
    [beta(ii),formresults] = run_reli(Prob_ii, Options);
    %beta(ii)               = run_reli_cMC(Prob_ii, Options);
    alphas{ii}             = formresults.alpha;
    flag_form(ii)          = formresults.flag;
end

% Commmented out because we do not have non-converged analyses and 
% `run_reli_subsetsim` is not harmonized with the rest of the code.
% for non-converged FORM analyses, do subset simulation
% nb: the reason FORM do not converge is because it is jumping between two
%     points as a result of the v_min check in EC2 shear formula
% if sum(flag_form) < n_ds
%     nonconvFORM = find(~flag_form);
%     n_ds_FORM   = length(nonconvFORM);
% %     beta_FORM   = beta(nonconvFORM);
%     print(['Number of non-converged FORM analyses: ', num2str(n_ds_FORM)])
%     for ii = 1:n_ds_FORM
%     %parfor ii = 1:n_ds_FORM
%         ds_num                                  = nonconvFORM(ii);
%         Prob_ii                                 = Prob(ds_num);
%         [beta(ds_num), ~] = run_reli_subsetsim(Prob_ii, Options);
%         alphas{ds_num}                          = {};
%     end
% %     beta_SS     = beta(nonconvFORM);
% %     beta_FORM - beta_SS
% end


% -------------------------------------------------------------------------
% Objective function
% -------------------------------------------------------------------------
%     fval   = trapz(chi_vec, weight.*(beta_t - Beta).^2);

weights = DS.weights_combis;
fval  = trapz(weights.*(beta_t - beta).^2);

% -------------------------------------------------------------------------
% Collect results
% -------------------------------------------------------------------------
if verbose == 0
   Results.alphas = alphas;
end
Results.beta            = beta;
Results.flag_form       = flag_form;

end