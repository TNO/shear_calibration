% plot Beta values versus chi1 and chi2 for given set of structural
% paramaters (3D surface plot)
%
function plot_3Dbeta(Results, DS, parset_of_interest)

% parameter set of interest
d_int    = parset_of_interest.d;
fcck_int = parset_of_interest.fcck;
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
combis_of_interest = combis(combis_set,[1,2,5]);
set_of_interest    = [d_int, fcck_int, rho_int];
logical_interest   = ismember(combis_of_interest,set_of_interest);
sum_matches        = sum(logical_interest,2);
sets_allmatches    = combis_set(sum_matches==3);
betas_of_interest  = Results.beta(sets_allmatches);
chi1_plot          = combis(sets_allmatches,3);
chi2_plot          = combis(sets_allmatches,4);

nchi1 = length(chi1_ds);
nchi2 = length(chi2_ds);
CHI1  = reshape(chi1_plot, nchi1, nchi2)';
CHI2  = reshape(chi2_plot, nchi1, nchi2)';
BETA  = reshape(betas_of_interest, nchi1, nchi2)';


% plot figure
scrsz = get(groot,'ScreenSize');  %2014version 
hFig = figure;
set(gcf,'PaperPositionMode','auto')
set(hFig, 'Position', [scrsz(3)/5 scrsz(3)/5 scrsz(3)/1.75 scrsz(4)/2.75])  % [left bottom width height]
set(gca,'units','normalized','position',[0.1 0.2 0.8 0.75]); 

subplot(1,2,1)
surf(CHI1,CHI2,BETA)
%alpha 0.5
%shading interp
%xlim([minX,maxX]);
%ylim([minY,maxY]);
%zlim([minZ,maxZ]);
%caxis([minZ,maxZ])
%colorbar  

%legend(legend_rvs, 'Location', 'SouthWest')
%title(['Estimated \alpha^{2} for: d = ' num2str(d_int) ', f_{cck} = ' ...
%       num2str(fcck_int) ', \rho_{s} = ' num2str(rho_int) ' and \chi_{2} = ' num2str(chi2_int)] )
xlabel('\chi_{1}')
ylabel('\chi_{2}')
zlabel('\beta')


subplot(1,2,2)
contourf(CHI1,CHI2,BETA)
xlabel('\chi_{1}')
ylabel('\chi_{2}')
h = colorbar;
set(get(h,'title'),'string','\beta');

