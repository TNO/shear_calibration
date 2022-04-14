clearvars; close all; clc

to_path = {'..\calibration\', '..\..\..\shared_utils\visualization\'};
cellfun(@(x) addpath(genpath(x)), to_path)

save_fig = 1;
beta_target = 0;
n = 100;
beta = beta_target + linspace(-0.5, 0.5, n);

fd1 = fdist_squared(beta, beta_target);
fd2 = fdist_hansen_sorensen(beta, beta_target);

% scale for better comparison
scaler = fd2(end) / fd1(end);

plot(beta - beta_target, scaler * fd1, 'black', 'LineWidth', 2)
hold on
plot(beta - beta_target, fd2, 'black--', 'LineWidth', 2)

xlabel('$\beta - \beta_\mathrm{target}$')
ylabel('$f_\mathrm{dist}$')
legend('squared', 'Hansen \& Sorensen')

prettify(gcf)

% SAVE
if save_fig == 1
    fwidth  = 10;
    fheight = 8;
    fpath   = '..\..\results\fdist_comparison';
    figuresize(fwidth , fheight , 'cm')
    export_fig(fpath, '-png', '-m2.5')
end

% Remove from the path: for safety as the functions added to the path are
% available from everywhere
cellfun(@(x) rmpath(genpath(x)), to_path)