% Testing the visualization

clearvars
% close all
clc

fname = 'tmp_ec2_weights';

load([fname, '.mat'])

% plot in Matlab
visu_calibr(Options, Results, DS)



