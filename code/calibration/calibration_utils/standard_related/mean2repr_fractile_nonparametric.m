% Calculate representative value from the meanand probabilistic and
% standardized models where the representative value is a fractile.
%
%SYNOPSYS
% [X, x_repr] = MEAN2REPR_FRACTILE_PARAMETRIC(x_mean, X)
%
%INPUT
% x_mean     Mean value of the distribution defined in `X` /scalar/.
% X          A struct
%
%OUTPUT
% [X, x_repr]
%
%NOTES
% Does not accept matrix inputs!
%


function [x_repr, X] = mean2repr_fractile_nonparametric(x_mean, X)

% -------------------------------------------------------------------------
% INITIALIZE
% -------------------------------------------------------------------------

P_repr = X.P_repr;
dist = X.dist;
dist_ID = X.dist_ID;

X.mean = x_mean;

% -------------------------------------------------------------------------
% CALCULATE REPR
% -------------------------------------------------------------------------

switch dist
    case 33
        [~, x_grid, cdf]    = hardcoded_cdf(1, dist_ID);
        x_repr = interp1(cdf, x_grid, P_repr);
    otherwise
        error(['Unkown or not implemented distributiont type (dist): ', num2str(dist)])
end

X.repr = x_repr;
X.repr_base = X.repr;

end