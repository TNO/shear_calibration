clearvars; close all; clc

theta_R = 1.0;
gamma_R = 1.0;

fc = linspace(20, 100, 500);
b = 100;
d = 300;
rho = 0.01;
d_lower = 16;
a_to_d_ratio = 3;

Asl = b * d * rho;

VR = MC2010_level_II_codified_2019(fc, Asl, b, d, d_lower, a_to_d_ratio, theta_R, gamma_R);

plot(fc, VR, 'LineWidth', 2)
xlabel('$f_\mathrm{c}$ [MPa]')
ylabel('$V_\mathrm{R}$ [kN]')
grid on

prettify(gcf)