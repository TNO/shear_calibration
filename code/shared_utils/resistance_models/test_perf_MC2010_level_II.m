% Test the vectorized implementation of the MC2010 level II resistance
% model + approximate timing
clearvars; close all; clc

theta_R = 1.0;
gamma_R = 1.0;

len_vecs = round(logspace(0, 4, 10));

% to activate the epsx > 0.003 part of the code as well
rhos = [0.01, 0.001];
b = 100;
d = 300;
d_lower = 16;
a_to_d_ratio = 3;

num_lens = length(len_vecs);

t_fzero = nan(1, num_lens);
t_vect = nan(1, num_lens);
for ii = 1:num_lens
    len_vec = len_vecs(ii);
    fc = linspace(20, 100, len_vec);
    disp('------------------------------------')
    disp(['vector length: ', num2str(len_vec)])

    for jj = 1:length(rhos)
        rho = rhos(jj);
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

        t_fzero(ii) = t_elapsed;
        t_vect(ii) = t_elapsed_vect;
        
        if ii == num_lens
            subplot(1,2,jj)
            plot(fc, VR, 'LineWidth', 2)
            hold on
            plot(fc, pVR, '--', 'LineWidth', 2)
            xlabel('$f_\mathrm{c}$ [MPa]')
            ylabel('$V_\mathrm{R}$ [kN]')
            grid on
            legend('loop \& fzero', 'vect', 'Location', 'southeast')
            title(['$\rho=$', sprintf('%.3e', rho)])
            
            prettify(gcf)
        end
    end
end

% compare runtime
figure
loglog(len_vecs, t_fzero, 'LineWidth', 2)
hold on
plot(len_vecs, t_vect, '--', 'LineWidth', 2)
legend('loop \& fzero', 'vect')
xlabel('Length of input vectors')
ylabel('Runtime [s]')
grid on
prettify(gcf)