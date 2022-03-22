% Calculate the mean from representative value when they are related 
% like this:
%
% x_mean = x_repr + shift

function [x_repr, X] = mean2repr_shift(x_mean, X, shift)

x_repr = x_mean - shift;
X.mean = x_mean;
X.repr = x_repr;

end