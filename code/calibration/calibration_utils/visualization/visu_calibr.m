% Visualize calibration results
%
%
% gramm is needed: https://nl.mathworks.com/matlabcentral/fileexchange/54465-gramm--complete-data-visualization-toolbox--ggplot2-r-like-

function visu_calibr(Options, Results, DS, fname, R_plot)

if nargin < 3
    fname = [];
end

if nargin < 4
    R_plot = 0;
end

d_ds            = DS.Range.d;
f_cck_ds        = DS.Range.f_cck;
chi1_ds         = DS.Range.chi1;
chi2_ds         = DS.Range.chi2;
rho_ds          = DS.Range.rho;

combis          = DS.combis;
combis_colnames = DS.combis_colnames;

% % create all possible combinations for the design scenarios
% combis          = combvec(d_ds, f_cck_ds, chi1_ds, chi2_ds, rho_ds)';
% combis_colnames = {'d', 'f_cck', 'chi1', 'chi2', 'rho'};
% 
% % calculate weights (based on chi1 and chi2):
% if isfield(Options, 'weights_filename') 
%     weights_chi1 = DS.weights_chi1;
%     weights_chi2 = DS.weights_chi2;
%     
%     chi1chi2 = combis(:,3:4);
%     for i = 1:length(chi1_ds)
%         chi1chi2(chi1chi2(:,1) == chi1_ds(i), 1) = weights_chi1(i);
%     end
%     for i = 1:length(chi2_ds)
%         chi1chi2(chi1chi2(:,2) == chi2_ds(i), 2) = weights_chi2(i);
%     end
%     DS.weights_combis = chi1chi2(:,1).*chi1chi2(:,2);
%     DF.weight_logic   = DS.weights_combis > 0;
% %     DF.weight_logic   = double(DS.weights_combis > 0);
% end

% mimic dataframe object
DF.beta       = Results.beta;
DF.d          = combis(:,1);
DF.f_cck      = combis(:,2);
DF.chi1       = combis(:,3);
DF.chi2       = combis(:,4);
DF.rho        = combis(:,5);
DF.load_comb  = DS.load_combs_all;


% mimic ggplot2
g = gramm('x',DF.chi1, 'y',DF.beta, 'color',DF.chi2, 'linestyle',DF.d);
% g = gramm('x',DF.chi1, 'y',DF.beta, 'color',DF.chi2, 'linestyle',DF.d, 'marker',DF.weight_logic);
g.facet_grid(DF.rho, DF.load_comb);
%g.facet_grid(DF.f_cck, DF.rho);
g.geom_point();
g.geom_line();
%g.set_names('row','$f_\mathrm{cc,k}$','column','$\rho_\mathrm{sl}$',...
g.set_names('row','$\rho_\mathrm{sl}$','column','$var. actions$',...
    'x','$\chi_1$','y','$\beta$',...
    'color','$\chi_2$','linestyle','$d$');
% g.set_names('row','$f_\mathrm{cc,k}$','column','$\rho_\mathrm{sl}$',...
%     'x','$\chi_1$','y','$\beta$',...
%     'color','$\chi_2$','linestyle','$d$','marker','weight');

figure('Position', [100, 100, 600*length(rho_ds), 200*length(f_cck_ds)]);
g.draw();

if exist('prettify', 'file') == 2
    prettify(gcf)
end

% DS.combis           = combis;
% DS.combis_colnames  = combis_colnames;

% -------------------------------------------------------------------------
% R plot
% -------------------------------------------------------------------------

% % if R_plot == 1
% %     % save the results to be used by R
% %     save_results_to_visu(fname, Options, Results, DS)
% % end

end