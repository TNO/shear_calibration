function int = weighted_integral(beta, DS, Options)

% .........................................................................
% Initialize
% .........................................................................
load_combs_all = DS.load_combs_all;
% p_ds_all = DS.p_ds_all;
p_ds_all_colnames = DS.p_ds_all_colnames;

p_ds_all_colnames_expected = {
    'chi1', 'chi2', 'd', 'f_cck', 'rho', 'd_lower', 'a_to_d_ratio'
};
assert(isequal(p_ds_all_colnames, p_ds_all_colnames_expected))

chi1 = DS.Range.chi1;
chi2 = DS.Range.chi2;

load_combs =  Options.load_combs;
beta_target = Options.beta_target;
objective_function = Options.objective_function;

switch lower(objective_function)
    case 'squared'
        fdist = @(x) fdist_squared(x, beta_target);
    case 'hansen_sorensen'
        fdist = @(x) fdist_hansen_sorensen(x, beta_target);
    otherwise
        error(['Unknown objective function:', resistance_model])
end

n_load_comb = length(load_combs);
n_p_ds = length(p_ds_all_colnames);

% .........................................................................
% Pre-process
% .........................................................................
n_discr_p_ds = nan(1, n_p_ds);
for ii = 1:n_p_ds
    n_discr_p_ds(ii) = length(DS.Range.(p_ds_all_colnames{ii}));
end

beta_mxs = cell(1, n_load_comb);
ints = nan(1, n_load_comb);
for ii = 1:n_load_comb
    load_comb_ii = load_combs{ii};
    chi_weights_ii = DS.Weight.(load_comb_ii).W';
    
    bm_ii = strcmp(load_combs_all, load_comb_ii);

    if strcmp(load_comb_ii, 'traffic')
        n_discr_p_ds_ii = n_discr_p_ds([1, 3:end]);
%         n_discr_nonchi_p_ds_ii = n_discr_p_ds_ii(2:end);
    else
        n_discr_p_ds_ii = n_discr_p_ds;
%         n_discr_nonchi_p_ds_ii = n_discr_p_ds_ii(3:end);
    end

    size_ii = num2cell(n_discr_p_ds_ii);
    beta_mx_ii = reshape(beta(bm_ii), size_ii{:});
    beta_mxs{ii} = beta_mx_ii;

    % .....................................................................
    % Numerical integral
    % .....................................................................
    % this automatically broadcast over the non-chi dimensions
    chi_weighted_beta_dist_mx_ii = chi_weights_ii .* fdist(beta_mx_ii);

    if strcmp(load_comb_ii, 'traffic')
        int_over_chi = trapz(chi1, chi_weighted_beta_dist_mx_ii, 1);
    else
        int_over_chi = trapz(chi2, trapz(chi1, chi_weighted_beta_dist_mx_ii, 1));
    end

    int_over_chi = squeeze(int_over_chi);

    % integrate over the non-chi dimensions - assuming that they all have
    % the same weight and they are 'categorical'
%     nonchi_weights = prod(1 ./ n_discr_nonchi_p_ds_ii);
%     int_ii = sum(int_over_chi(:) * nonchi_weights);
    int_ii = sum(int_over_chi(:));
    ints(ii) = int_ii;
end

int = sum(ints);

end