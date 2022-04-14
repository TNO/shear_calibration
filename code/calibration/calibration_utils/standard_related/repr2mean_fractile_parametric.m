% Calculate mean from representative value and probabilistic and
% standardized models where the representative value is a fractile.
%
%SYNOPSYS
% [x_mean, X] = REPR2MEAN_FRACTILE_PARAMETRIC(x_repr, X)
%
%INPUT
% x_repr     Representative value that is a fractile of the distribution
%            defined in `X` /scalar/.
% X          A struct that must have at least the following fields:
%   .P_repr    P_repr = F(x_repr) /scalar/.
%   .cov       Coefficient of variation /scalar/. One and only of `.cov`
%              and `.std` must be NaN.
%   .std       Standard deviation /scalar/.
%   .dist      FERUM like distribution type indicator /integer/.
%
%OUTPUT
% [x_mean, X]
%
%NOTES
% Does not accept matrix inputs!
% WARNING: assumes mean of 1 for the base function of par-nonpar
% distributions!
%
%
% TODO:
% overkill


function [x_mean, X] = repr2mean_fractile_parametric(x_repr, X)

% -------------------------------------------------------------------------
% INITIALIZE
% -------------------------------------------------------------------------

P_repr      = X.P_repr;
cov         = X.cov;
std         = X.std;
dist        = X.dist;

if (~isnan(cov) && ~isnan(std)) || (isnan(cov) && isnan(std))
    error('Only one of the `X.cov` `X.std` pairs must differ from NaN.')
end

X.repr = x_repr;

% -------------------------------------------------------------------------
% CALCULATE MEAN
% -------------------------------------------------------------------------

switch dist
    case 1
        if ~isnan(cov)
            x_mean = x_repr/(norminv(P_repr)*cov + 1);
        else
            x_mean = x_repr - std * norminv(P_repr);
        end
        X.mean = x_mean;
    case 2
        if ~isnan(cov)
            x_mean = fzero(@(t) my_lognormcdf(x_repr, t, cov) - P_repr, x_repr);
        else
            x_mean = fzero(@(t) my_lognormcdf(x_repr, t, std./t) - P_repr, x_repr);
        end
        X.mean = x_mean;
    case 11
        if ~isnan(cov)
            x_mean = fzero(@(t) my_gumbelcdf(x_repr, t, cov) - P_repr, x_repr);
        else
            x_mean = fzero(@(t) my_gumbelcdf(x_repr, t, std./t) - P_repr, x_repr);
        end
        X.mean = x_mean;
    otherwise
        error(['Unkown or not implemented distributiont type (dist): ', num2str(dist)])
end

end