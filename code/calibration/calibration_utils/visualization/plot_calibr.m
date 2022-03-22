function plot_calibr(alpha_R, Prob, Options)

save_plot = 1;

%--------------------------------------------------------------------------
% PRE-PROCESSING
%--------------------------------------------------------------------------
chi_vec      = Options.chi_vec;
beta_t       = Options.beta_t;
alpha_R_ecov = Options.ecov.alpha_R;
ID           = Options.ID;
n_chi        = length(chi_vec);
cmp          = lines(8);

chi_vec      = chi_vec(:);


% get the beta value
% overkill, but it is still cheap
[~, Results_0]  = obj_fun(alpha_R_ecov, Prob, Options);
[~, Results_c]  = obj_fun(alpha_R, Prob, Options);


%--------------------------------------------------------------------------
% PLOT
%--------------------------------------------------------------------------

%..........................................................................
% RELIABILITY LEVEL
%..........................................................................
beta_0      = Results_0.beta;
beta_c      = Results_c.beta;

figure('Position', [100, 200, 700, 300])
plot(chi_vec, beta_0, 'Color', cmp(1,:))
hold on
plot(chi_vec, beta_c, 'Color', cmp(2,:))

plot(chi_vec, beta_t*ones(1,n_chi), '--red')

xx = min(xlim) + 0.8*range(xlim);
yy = beta_0(end) + 0.03*range(ylim);
text(xx, yy, ['$\alpha_\mathrm{R}=', num2str(round(alpha_R_ecov,3)), '$'], 'Color', cmp(1,:))

xx = min(xlim) + 0.8*range(xlim);
yy = beta_c(end) + 0.03*range(ylim);
text(xx, yy, ['$\alpha_\mathrm{R}=', num2str(round(alpha_R,3)), '$'], 'Color', cmp(2,:))

xx = min(xlim) + 0.05*range(xlim);
yy = beta_t + 0.03*range(ylim);
text(xx, yy, ['$\beta_\mathrm{target}=', num2str(round(beta_t,3)), '$'], 'Color', 'red')

xlabel('$\chi = Q_\mathrm{k}/(G_\mathrm{k} + Q_\mathrm{k})$')
ylabel('$\beta$')
leg = legend('original E-COV','calibrated E-COV');
leg.Location = 'eastoutside';
prettify(gcf)

if save_plot == 1
    set(gcf, 'Color', 'w');
    set(gca,'Fontsize', 10);
    figuresize(16 , 6 , 'cm')
    fpath   = [pwd,'\figures\',ID,'_beta_chi.png'];
    %     hgexport(gcf, fpath, s);
    export_fig(fpath, '-png', '-m2.5')
end

%..........................................................................
% DESIGN VALUE
%..........................................................................
R_d_ecov_0      = Results_0.R_d_ecov*ones(n_chi,1);
R_d_reli_0      = Results_0.R_d_reli;
R_d_ecov_c      = Results_c.R_d_ecov*ones(n_chi,1);
R_d_reli_c      = Results_c.R_d_reli;

figure('Position', [100, 200, 700, 300])
plot(chi_vec, R_d_ecov_0, '-', 'Color', cmp(1,:))
hold on
plot(chi_vec, R_d_reli_0, '--', 'Color', cmp(1,:))

plot(chi_vec, R_d_ecov_c, '-', 'Color', cmp(2,:))
plot(chi_vec, R_d_reli_c, '--', 'Color', cmp(2,:))

xx = min(xlim) + 0.8*range(xlim);
yy = R_d_ecov_0(end) + 0.03*range(ylim);
text(xx, yy, ['$\alpha_\mathrm{R}=', num2str(round(alpha_R_ecov,3)), '$'], 'Color', cmp(1,:))

xx = min(xlim) + 0.8*range(xlim);
yy = R_d_ecov_c(end) + 0.03*range(ylim);
text(xx, yy, ['$\alpha_\mathrm{R}=', num2str(round(alpha_R,3)), '$'], 'Color', cmp(2,:))


xlabel('$\chi = Q_\mathrm{k}/(G_\mathrm{k} + Q_\mathrm{k})$')
ylabel('$R_\mathrm{d}$')
leg = legend('original E-COV', 'reli. analysis (original E-COV)', 'calibrated E-COV', 'reli. analysis (calibrated E-COV)');
leg.Location = 'eastoutside';
prettify(gcf)

if save_plot == 1
    set(gcf, 'Color', 'w');
    set(gca,'Fontsize', 10);
    figuresize(16 , 6 , 'cm')
    fpath   = [pwd,'\figures\',ID,'_R_d_chi.png'];
    %     hgexport(gcf, fpath, s);
    export_fig(fpath, '-png', '-m2.5')
end

%..........................................................................
% ALPHA_R
%..........................................................................
alpha_R_ecov_0      = Options.ecov.alpha_R*ones(n_chi,1);
alpha_R_reli_0      = Results_0.alpha_R_reli;
alpha_R_ecov_c      = alpha_R*ones(n_chi,1);
alpha_R_reli_c      = Results_c.alpha_R_reli;

figure('Position', [100, 200, 700, 300])
plot(chi_vec, alpha_R_ecov_0, '-', 'Color', cmp(1,:))
hold on
plot(chi_vec, alpha_R_reli_0, '--', 'Color', cmp(1,:))

plot(chi_vec, alpha_R_ecov_c, '-', 'Color', cmp(2,:))
plot(chi_vec, alpha_R_reli_c, '--', 'Color', cmp(2,:))

xx = min(xlim) + 0.8*range(xlim);
yy = alpha_R_ecov_0(end) + 0.03*range(ylim);
text(xx, yy, ['$\alpha_\mathrm{R}=', num2str(round(alpha_R_ecov,3)), '$'], 'Color', cmp(1,:))

xx = min(xlim) + 0.8*range(xlim);
yy = alpha_R_ecov_c(end) + 0.03*range(ylim);
text(xx, yy, ['$\alpha_\mathrm{R}=', num2str(round(alpha_R,3)), '$'], 'Color', cmp(2,:))

ylim([0,1])
xlabel('$\chi = Q_\mathrm{k}/(G_\mathrm{k} + Q_\mathrm{k})$')
ylabel('$\alpha_{R}$')
leg = legend('original E-COV', 'reli. analysis (original E-COV)', 'calibrated E-COV', 'reli. analysis (calibrated E-COV)');
leg.Location = 'eastoutside';
prettify(gcf)

if save_plot == 1
    set(gcf, 'Color', 'w');
    set(gca,'Fontsize', 10);
    figuresize(16 , 6 , 'cm')
    fpath   = [pwd,'\figures\',ID,'_alpha_R_chi.png'];
    %     hgexport(gcf, fpath, s);
    export_fig(fpath, '-png', '-m2.5')
end

end
