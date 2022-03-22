% Estimate the model uncertainty (C_c) of a selected shear resistance formula.
%
% kappa = V_Rexp./V_Rmodel
% kappa is assuemed to be lognormally distributed
% the parameters of the lognormal distribution are estimated using the maximum
% likelihood method
%
%
% TODO: check and fix model = 'mc2010-current'; probably there is a problem
% with units.

clearvars
close all
clc

% load functions from shared utilities (Matlab...)
to_path = '..\..\shared_utils\resistance_models\';
addpath(genpath(to_path))

% =========================================================================
% OPTIONS
% =========================================================================

model       = 'ec2-current';
% model       = 'ec2-proposed';
% model       = 'mc2010-current';

consider_vrmin = false;
% confidence level
ci          = 0.95;

main_title  = ['Calibrated ', model];
sizedata    = 30;

fmat        = 'beam_p_adj.mat';

save_fig    = 1;

% =========================================================================
% LOAD AND PREPROCESS DATA
% =========================================================================

fpath       = ['..\..\..\data\', fmat];
load(fpath);

V_Rexp      = Vu;
gamma_C     = 1.0;         % the partial safety factor for concrete is 1.5 (for our purpose, it is excluded!)
fc          = fc;          % in [MPa], according to par. 7.2.3 (for our purpose, we use the fc;test value; otherwise: fck = fcmean - 8)
b           = bw;          % in [mm], width of the beam
d           = d;           % in [mm], effective height
rho         = rho/100;
dg          = 16*ones(size(rho)); % note that this done as dg is fixed in the calibration!!
Asl         = rho.*b.*d;          % in [mm2], area of tensile reinforcement in considered section

% =========================================================================
% ESTIMATE parameters
% =========================================================================
% multiplicative error
switch lower(model)
    case 'ec2-current'
        shear_formula     = @(fc, Asl, b, d, C_c, gamma_C) EC2_codified_2019(fc, Asl, b, d, C_c, gamma_C, consider_vrmin);
    case 'ec2-proposed'
        shear_formula     = @(fc, Asl, b, d, C_c, gamma_C) EC2_proposed_Yuguang_2019(fc, Asl, b, d, C_c, gamma_C);
%         shear_formula     = @(fc, Asl, b, d, C_c, gamma_C) EC2_proposed_TG4_2016(fc, Asl, b, d, C_c, gamma_C);
    case 'mc2010-current'
        shear_formula     = @(fc, Asl, b, d, C_c, gamma_C) MC2010_level_II_codified_2019(fc, Asl, b, d, C_c, gamma_C);
    otherwise
        error('Unknown model.')
end

V_Rmodel                = shear_formula(fc, Asl, b, d, 1, 1);
kappa                   = V_Rexp./V_Rmodel;

% simplified 
[cpar, min_nLL]         = fit_lognorm2_mle(kappa, 'par');
[c_mean, c_cov]         = lognormstat(cpar(1), cpar(2), 'par');

C_c                     = c_mean;
V_Rmod                  = shear_formula(fc, Asl, b, d, C_c, gamma_C);

disp('Model uncertainty maximum likelihood estimate, C_c')
disp(['     mean of C_c : ', sprintf('%.3f', c_mean)])
disp(['     cov of C_c  : ', sprintf('%.3f', c_cov)])
fprintf('\n')

% -------------------------------------------------------------------------
% Goodness-of-fit measures
% -------------------------------------------------------------------------

np                      = length(V_Rmod);

% mean absolute error
MAE                     = mean(abs(V_Rexp - V_Rmod));

% median absolute error
MEDAE                   = median(abs(V_Rexp - V_Rmod));

% AIC
AIC                     = 2*1 + 2*min_nLL;

% RMSD
RMSD                    = sqrt(sum((V_Rexp - V_Rmod).^2)/np);

% rho
rho_c                   = corr(V_Rexp, V_Rmod);

GoF                     = table(MAE, MEDAE, AIC, RMSD, rho_c);

disp(GoF)

% =========================================================================
% VISUALIZE
% =========================================================================

% -------------------------------------------------------------------------
% Ratio Vexp/VRmodcalibr
% -------------------------------------------------------------------------

figure

hs                      = scatter(d, V_Rexp./V_Rmod);
hs.MarkerEdgeColor      = 1*ones(1,3);
hs.MarkerFaceColor      = 0.5*ones(1,3);
hs.MarkerFaceAlpha      = 0.5;
hs.SizeData             = sizedata;

hold on

xx = linspace(min(d), max(d), 1e2);
plot(xx, ones(size(xx)), 'black')

xlabel('Effective depth, $d$ [mm]')
ylabel('$V_\mathrm{R,exp}/V_\mathrm{R,model,mean}$')
title(main_title)
prettify(gcf)

% SAVE
if save_fig == 1
    fwidth  = 10;
    fheight = 10;
    fpath   = ['./results/',model,'_exp_mod_ratio'];
    figuresize(fwidth , fheight , 'cm')
    export_fig(fpath, '-png', '-m2.5')
end

% -------------------------------------------------------------------------
% Observed vs Predicted
% -------------------------------------------------------------------------
figure
hs                  = scatter(V_Rmod, V_Rexp);
hs.MarkerEdgeColor  = 1*ones(1,3);
hs.MarkerFaceColor  = 0.5*ones(1,3);
hs.MarkerFaceAlpha  = 0.5;
hs.SizeData         = sizedata;
hold on

VV          = [V_Rmod(:); V_Rexp(:)];
xx          = linspace(min(VV), max(VV), 1e2);
plot(xx, xx, 'black')
mm          = C_c;
mu          = lognorminv((1+ci)/2, C_c, c_cov)/C_c;
ml          = lognorminv((1-ci)/2, C_c, c_cov)/C_c;

plot(xx, xx*mu, 'black--')
plot(xx, xx*ml, 'black--')

axis equal
xlim([0, max(mm(:), max(VV))])
ylim([0, max(mm(:), max(VV))])


xlabel('Predicted, $V_\mathrm{R,model,mean}$ [kN]')
ylabel('Observed, $V_\mathrm{R,exp}$ [kN]')
title([main_title, '; confidence level: ', num2str(ci)])
prettify(gcf)

% SAVE
if save_fig == 1
    fwidth  = 10;
    fheight = 10;
    fpath   = ['./results/',model,'_exp_vs_pred'];
    figuresize(fwidth , fheight , 'cm')
    export_fig(fpath, '-png', '-m2.5')
end

% -------------------------------------------------------------------------
% Residual
% -------------------------------------------------------------------------
figure('Position', [100, 100, 800, 400])
subplot(1,2,1)
hs                  = scatter(1:np, V_Rmod-V_Rexp);
hs.MarkerEdgeColor  = 1*ones(1,3);
hs.MarkerFaceColor  = 0.5*ones(1,3);
hs.MarkerFaceAlpha  = 0.5;
hs.SizeData         = sizedata;
hold on

plot(1:np, zeros(size(V_Rmod)), 'black')
title(main_title)
xlabel('Experiment')
ylabel('Residual error, $V_\mathrm{R,model,mean} - V_\mathrm{R,exp}$ [kN]')
ylim([-50, 150])
prettify(gcf)

% -------------------------------------------------------------------------
% Histogram
% -------------------------------------------------------------------------
subplot(1,2,2)
cc          = V_Rexp./shear_formula(fc, Asl, b, d, 1, gamma_C);

[f, x]      = hist(cc);
f           = f./trapz(x,f);

hb                  = bar(x,f);
hb.BarWidth         = 1;
hb.FaceColor        = 0.5*ones(1,3);
hb.FaceAlpha        = 0.5;
hb.EdgeColor        = 1*ones(1,3);

hold on
xx = linspace(min(cc), max(cc), 1e2);
yy = lognormpdf(xx, c_mean, c_cov);
plot(xx, yy, 'black')

box off
set(gca,'YColor',[1 1 1]);
xlabel('$C$ [-]')
title(main_title)
prettify(gcf)

% SAVE
if save_fig == 1
    fwidth = 18;
    fheight = 8;
    fpath   = ['./results/',model,'_resi_and_C_histogram'];
    figuresize(fwidth , fheight , 'cm')
    export_fig(fpath, '-png', '-m2.5')
end

% Remove from the path: for safety as the functions added to the path are
% available from everywhere
rmpath(genpath(to_path))