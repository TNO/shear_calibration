%TODO
% * collect all formresults into an array of structure
%
%

function [fval, Results, DS] = obj_fun(x, Prob, Prob_actions, DS, Options)

% -------------------------------------------------------------------------
% Initialize
% -------------------------------------------------------------------------
verbose                 = Options.verbose;

% -------------------------------------------------------------------------
% Semi-probabilistic design
% -------------------------------------------------------------------------
% semi-probabilistic design
free_par                = x;
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
    [beta(ii), formresults] = run_reli(Prob_ii, Options, 'form');
    
    alphas{ii}             = formresults.alpha;
    converged(ii)          = formresults.converged;

    if formresults.converged == 0
        [beta(ii), simulationresults] = run_reli(Prob_ii, Options, 'is', formresults);
        
        if simulationresults.converged == 1 && verbose > 1
            disp('Importance sampling has converged with the last step of not converged FORM.')
        end
        alphas{ii}             = nan(size(formresults.alpha));
        converged(ii)          = simulationresults.converged;
    end
    
end

% num_no_conv_form_is = sum(~converged);
% For non-converged FORM + IS analyses we perform subset simulation
% nb: FORM does not converge because it is jumping between 
% two points as a result of the v_min check in the EC2 shear formula
% WATCHOUT: we do not have a convergence check on SS (we should not reach
% SS)
% if num_no_conv_form_is > 0
%     disp(['Number of not converged FORM + IS analyses per objective function ' ...
%         'evaluation: ', num2str(num_no_conv_form_is)])
%     disp('Subset simulation is being run for these not converged analyses...')
%     no_conv_form_is_idx = find(~converged);
%     beta_no_conv_form_is   = beta(no_conv_form_is_idx);
%     
%     beta_ss = nan(num_no_conv_form_is, 1);
%     for ii = 1:num_no_conv_form_is
% %     parfor ii = 1:num_no_conv_form
%         ds_idx  = no_conv_form_is_idx(ii);
%         Prob_ii = Prob(ds_idx);
%         [beta_ss(ii), ~] = run_reli(Prob_ii, Options, 'subset');
%     end
% 
%     % replace the non-converged FORM results with SS results
%     beta(no_conv_form_is_idx) = beta_ss;
%     beta_summary_table = table(beta_no_conv_form_is, beta_ss);
%     disp(beta_summary_table)
% end


% -------------------------------------------------------------------------
% Objective function
% -------------------------------------------------------------------------

% beta_t = Options.beta_target;
% weights = DS.weights_combis;
% fval  = trapz(weights.*(beta_t - beta).^2);
fval = weighted_integral(beta, DS, Options);

% -------------------------------------------------------------------------
% Collect results
% -------------------------------------------------------------------------

Results.alphas          = alphas;
Results.beta            = beta;
Results.converged       = converged;

end