% Save results for visualization (in R)
%
%
%

function save_results_for_visu(fname, output)

load(['./results/', fname, '.mat'], 'Results', 'DS')

% create all possible combinations for the design scenarios
p_ds_all        = DS.p_ds_all;

% weight_logic    = DS.weights_combis > 0;
weight          = DS.weights_combis;
load_comb       = DS.load_combs_all;

chi1    = p_ds_all(:,1);
chi2    = p_ds_all(:,2);
d       = p_ds_all(:,3);
f_cck   = p_ds_all(:,4);
rho     = p_ds_all(:,5);
d_lower = p_ds_all(:,6);
a_to_d_ratio = p_ds_all(:,7);

% mimic dataframe object
switch lower(output)
    case 'beta'
        beta    = Results.beta;
        % beta    = round(beta*1e4)/1e4;
        
        T       = table(beta, load_comb, d, f_cck, chi1, chi2, rho, d_lower, a_to_d_ratio, weight);
        
        tpath   = ['./results/beta_', fname, '.csv'];
        writetable(T, tpath)
        
    case 'alpha'
        rvs  = {'theta_R', 'f_c', 'd', 'b', 'A_sl', 'G', 'theta_G', 'Q1', 'theta_Q1', 'Q2', 'theta_Q2', 'theta_E'};
        cols = length(rvs);
        rows = length(Results.alphas);
        
        alpha = zeros(rows,cols);
        
        % to be optimized: alphas_i currently not same length (due to
        % traffic load combi --> no Q2). 
        for i = 1:rows
            alpha_i       = (Results.alphas{i})';
            n             = length(alpha_i); 
            alpha(i,1:n) = alpha_i;
        end
        
        alpha_r      = reshape(alpha,rows*cols,1);
        alpha2_r     = alpha_r.^2;
        alpha_labels = reshape(repmat(rvs,rows,1),rows*cols,1);       
        
        load_comb_r = repmat(load_comb,cols,1);
        d_r         = repmat(d,cols,1);
        f_cck_r     = repmat(f_cck,cols,1);
        chi1_r      = repmat(chi1,cols,1);
        chi2_r      = repmat(chi2,cols,1);
        rho_r       = repmat(rho,cols,1);

        d_lower_r   = repmat(d_lower, cols, 1);
        a_to_d_ratio_r   = repmat(a_to_d_ratio, cols, 1);
        weight_r    = repmat(weight,cols,1);

        T        = table(alpha2_r, alpha_r, alpha_labels, load_comb_r, d_r, f_cck_r, chi1_r, chi2_r, rho_r, d_lower_r, a_to_d_ratio_r, weight_r);
        
        tpath   = ['./results/alpha2_', fname, '.csv'];
        writetable(T, tpath)
        
    otherwise
        error(['Unknown output item: ', output])
end
        
disp('The results are sucessfully written to the hard drive.')

end