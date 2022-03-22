% Calculate representative value from the meanand probabilistic and
% standardized models where the representative value is a fractile.
%
%SYNOPSYS
% [X, x_repr] = MEAN2REPR_FRACTILE_PARAMETRIC(x_mean, X)
%
%INPUT
% x_mean     Mean value of the distribution defined in `X` /scalar/.
% X          A struct that must have at least the following fields:
%   .P_repr    P_repr = F(x_repr) /scalar/.
%   .cov       Coefficient of variation /scalar/. One and only of `.cov`
%              and `.std` must be NaN.
%   .std       Standard deviation /scalar/.
%   .dist      FERUM like distribution type indicator /integer/.
%
%OUTPUT
% [X, x_repr]
%
%NOTES
% Does not accept matrix inputs!
%
%
% TODO:
% overkill


function [x_repr, X] = mean2repr_fractile_parametric(x_mean, X)

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

if ~isnan(std)
    cov = std / x_mean;
end

% -------------------------------------------------------------------------
% CALCULATE REPR
% -------------------------------------------------------------------------

switch dist
    case 1
        x_repr = norminv(P_repr)*std + x_mean;
    case 2
        x_repr = fzero(@(t) my_lognormcdf(t, x_mean, cov) - P_repr, x_mean);
    case 11
        x_repr = fzero(@(t) my_gumbelcdf(t, x_mean, cov) - P_repr, x_mean);
    otherwise
        error(['Unkown or not implemented distributiont type (dist): ', num2str(dist)])
end

X.repr = x_repr;

end