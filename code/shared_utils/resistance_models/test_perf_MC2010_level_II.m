% 
clearvars; close all; clc

theta_R = 1.0;
gamma_R = 1.0;

fc = linspace(20, 100, 5000);
% to activate the epsx > 0.003 part of the code as well
rhos = [0.01, 0.001];
b = 100;
d = 300;
d_lower = 16;
a_to_d_ratio = 3;

for ii = 1:length(rhos)
    rho = rhos(ii);
    Asl = b * d * rho;
    
    t_start = tic;
    VR = fzero_MC2010_level_II_codified_2019( ...
        fc, Asl, b, d, d_lower, a_to_d_ratio, theta_R, gamma_R);
    t_elapsed = toc(t_start);
    
    t_start = tic;
    pVR = MC2010_level_II_codified_2019( ...
        fc, Asl, b, d, d_lower, a_to_d_ratio, theta_R, gamma_R);
    t_elapsed_vect = toc(t_start);
    
    assert(all(abs(pVR - VR) < 1e-4))
    
    disp(['root finding & loop (fzero):  ', sprintf('%.3e', t_elapsed), ' sec'])
    disp(['vectorized:                   ', sprintf('%.3e', t_elapsed_vect), ' sec'])
    
    subplot(1,2,ii)
    plot(fc, VR, 'LineWidth', 2)
    hold on
    plot(fc, pVR, '--', 'LineWidth', 2)
    xlabel('$f_\mathrm{c}$ [MPa]')
    ylabel('$V_\mathrm{R}$ [kN]')
    grid on
    title(['$\rho=$', sprintf('%.3e', rho)])
    
    prettify(gcf)
end