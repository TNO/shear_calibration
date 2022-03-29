% Save results for visualization (in R)
%
%
%

function save_results_for_visu(fname, output)

load(['./results/', fname, '.mat'], 'Results', 'DS')

% create all possible combinations for the design scenarios
combis          = DS.combis;

% weight_logic    = DS.weights_combis > 0;
weight          = DS.weights_combis;
load_comb       = DS.load_combs_all;

d       = combis(:,1);
f_cck   = combis(:,2);
chi1    = combis(:,3);
chi2    = combis(:,4);
rho     = combis(:,5);
d_lower = combis(:,6);
a_to_d_ratio = combis(:,7);

% mimic dataframe object
switch lower(output)
    case 'beta'
        beta    = Results.beta;
        % beta    = round(beta*1e4)/1e4;
        
        T       = table(beta, load_comb, d, f_cck, chi1, chi2, rho, d_lower, a_to_d_ratio, weight);
        
        tpath   = ['./results/beta_', fname, '.csv'];
        writetable(T, tpath)
        
    case 'alpha'
        rvs  = {'theta_R', 'f_cck', 'd', 'b', 'rho_s', 'G', 'theta_G', 'Q1', 'theta_Q1', 'Q2', 'theta_Q2', 'theta_E'};
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
       
%         alpha2_C   = alpha2(:,1);
%         alpha2_fcc = alpha2(:,2);
%         alpha2_d   = alpha2(:,3);
%         alpha2_b   = alpha2(:,4);
%         alpha2_rho = alpha2(:,5);
%         alpha2_G   = alpha2(:,6);
%         alpha2_KG  = alpha2(:,7);
%         alpha2_Q1  = alpha2(:,8);
%         alpha2_KQ1 = alpha2(:,9);
%         alpha2_Q2  = alpha2(:,10);
%         alpha2_KQ2 = alpha2(:,11);
        
%         T       = table(alpha2_C, alpha2_fcc, alpha2_d, alpha2_b, alpha2_rho, ...
%                         alpha2_G, alpha2_KG, alpha2_Q1, alpha2_KQ1, alpha2_Q2, ...
%                         alpha2_KQ2, load_comb, d, f_cck, chi1, chi2, rho, weight);
        
        tpath   = ['./results/alpha2_', fname, '.csv'];
        writetable(T, tpath)
        
    otherwise
        error(['Unknown output item: ', output])
end
        
disp('The results are sucessfully written to the hard drive.')

end