% Calculate the representative value from the mean when they are related 
% like this:
%
% x_mean = x_repr + shift

function [x_mean, X] = repr2mean_shift(x_repr, X, shift)

x_mean = x_repr + shift;
X.mean = x_mean;
X.repr = x_repr;

end