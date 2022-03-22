% plot alpha - values versus chi1 and chi2 for given set of structural
% paramaters
%
function plot_alphas(Results, DS, parset_of_interest)

% parameter set of interest
d_int    = parset_of_interest.d;
fcck_int = parset_of_interest.fcck;
chi1_int = parset_of_interest.chi1;
chi2_int = parset_of_interest.chi2;
rho_int  = parset_of_interest.rho;
set_int  = parset_of_interest.set;


% initialization
d_ds            = DS.Range.d;
f_cck_ds        = DS.Range.f_cck;
chi1_ds         = DS.Range.chi1;
chi2_ds         = DS.Range.chi2;
rho_ds          = DS.Range.rho;

legend_rvs      = {'C', 'f_{cck}', 'd', 'b', '\rho_{s}', 'G', 'KG', 'Q1', 'KQ1', 'Q2', 'KQ2'};

combis          = DS.combis;
combis_colnames = DS.combis_colnames;

% % create all possible combinations for the design scenarios
% combis          = combvec(d_ds, f_cck_ds, chi1_ds, chi2_ds, rho_ds)';
% combis_colnames = {'d', 'f_cck', 'chi1', 'chi2', 'rho'};

combis_set      = find(ismember(DS.load_combs_all, set_int));


% get alpha-values in combis for alpha - chi1 plot and alpha - chi2 plot
combis_chi1plot = round(combis(combis_set,[1:2,4:5]),4);
combi_int_chi1  = [d_int, fcck_int, chi2_int, rho_int];
logical_chi1    = ismember(combis_chi1plot,combi_int_chi1);
sum_matches     = sum(logical_chi1,2);
chi1_allmatches = combis_set(sum_matches==4);
for i = 1:length(chi1_allmatches)
    alphas2_chi1plot(i,:) = (Results.alphas{chi1_allmatches(i)}).^2;
end

combis_chi2plot = round(combis(combis_set,[1:3,5]),4);
combi_int_chi2  = [d_int, fcck_int, chi1_int, rho_int];
logical_chi2    = ismember(combis_chi2plot,combi_int_chi2);
sum_matches     = sum(logical_chi2,2);
chi2_allmatches = combis_set(sum_matches==4);
for i = 1:length(chi2_allmatches)
    alphas2_chi2plot(i,:) = (Results.alphas{chi2_allmatches(i)}).^2;
end


% plot figures
figure
area(chi1_ds', alphas2_chi1plot)
ylim([0,1])

legend(legend_rvs, 'Location', 'SouthWest')

title(['Estimated \alpha^{2} for: d = ' num2str(d_int) ', f_{cck} = ' ...
       num2str(fcck_int) ', \rho_{s} = ' num2str(rho_int) ' and \chi_{2} = ' num2str(chi2_int)] )
xlabel('\chi_{1}')
ylabel('\alpha^{2}')


figure
area(chi2_ds', alphas2_chi2plot)
ylim([0,1])

legend(legend_rvs, 'Location', 'SouthWest')

title(['Estimated \alpha^{2} for: d = ' num2str(d_int) ', f_{cck} = ' ...
       num2str(fcck_int) ', \rho_{s} = ' num2str(rho_int) ' and \chi_{1} = ' num2str(chi1_int)] )
xlabel('\chi_{2}')
ylabel('\alpha^{2}')
