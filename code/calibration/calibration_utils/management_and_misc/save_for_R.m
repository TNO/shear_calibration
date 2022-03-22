clearvars
close all
clc

% ID = '2018-Jul-04_13.53.55';
% ID = '2018-Jul-15_13.49.31';
% ID = '2018-Jul-15_16.12.29';
% ID = '2018-Jul-19_17.37.23';
ID = '2018-Jul-19_15.08.36';

output = {'beta', 'alpha'};

for i = 1:length(output)
    fname  = ['calibration_run_', ID];
    save_results_for_visu(fname, output{i})
end

% fname = 'calibration_run_2018-Jul-05_11.15.27';
% save_results_for_visu(fname)