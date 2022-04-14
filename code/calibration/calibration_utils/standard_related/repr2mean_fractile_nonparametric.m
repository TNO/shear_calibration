% Calculate mean from representative value and probabilistic and
% standardized models where the representative value is a fractile.
%
%SYNOPSYS
% [x_mean, X] = REPR2MEAN_FRACTILE_NONPARAMETRIC(x_repr, X)
%
%INPUT
% x_repr     Representative value that is a fractile of the distribution
%            defined in `X` /scalar/.
% X          A struct

%
%OUTPUT
% [x_mean, X]
%
%NOTES
% Does not accept matrix inputs!
% WARNING: assumes mean of 1 for the base function of par-nonpar
% distributions!


function [x_mean, X] = repr2mean_fractile_nonparametric(x_repr, X)

% -------------------------------------------------------------------------
% INITIALIZE
% -------------------------------------------------------------------------

dist = X.dist;

X.repr = x_repr;

% -------------------------------------------------------------------------
% CALCULATE MEAN
% -------------------------------------------------------------------------

switch dist
    case 33
        mean_base   = X.mean_base;
        repr_base   = X.repr_base;

        x_scale     = x_repr/repr_base;
        x_mean      = mean_base*x_scale;

        X.mean      = x_mean;
        X.scale     = x_scale;
        X.shift     = 0;
    otherwise
        error(['Unkown or not implemented distributiont type (dist): ', num2str(dist)])
end

end