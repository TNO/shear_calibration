% Calculate mean from a fractile (characteristic value typically)
%
%SYNOPSYS
% X = CHAR2MEAN(x_char, X)
%
%INPUT
% x_frac    fractile value /scalar/
% X
%   .P         P_frac = F(x_frac) /scalar/
%   .cov       coefficient of variation /scalar/
%   .dist      FERUM like distribution type indicator /integer/
%
%OUTPUT
% X
%
%NOTES
% Does not accept matrix inputs!
% WARNING: assumes mean of 1 for the base function of par-nonpar
% distributions!
%
%
% TODO:
% overkill


function [X, x_mean] = char2mean(x_char, X)

% -------------------------------------------------------------------------
% INITIALIZE
% -------------------------------------------------------------------------

P_char      = X.P_char;
cov         = X.cov;
dist        = X.dist;
dist_ID     = X.dist_ID;

X.char = x_char;

% -------------------------------------------------------------------------
% CALCULATION
% -------------------------------------------------------------------------

switch dist
    case 1
        x_mean      = x_char/(norminv(P_char)*cov + 1);
        X.mean      = x_mean;
    case 2
        x_mean      = fzero(@(t) my_lognormcdf(x_char, t, cov) - P_char, x_char);
        X.mean      = x_mean;
    case 11
        x_mean      = fzero(@(t) my_gumbelcdf(x_char, t, cov) - P_char, x_char);
        X.mean      = x_mean;
    case 32
        clear nonparametric_cdf
        
        mean_base   = X.mean_base;
        char_base   = X.char_base;
% 
%         x0          = [1e-4, 10];
%         f1          = my_parnonparfun(x0(1));
%         f2          = my_parnonparfun(x0(2));
%         while isnan(f1)
%             x0(1) = x0(1) + 0.01;
%             f1          = my_parnonparfun(x0(1));
%         end
%         
%         if sign(f1) == sign(f2)
%             pi
%         end
%         
%         x_scale     = fzero(@my_parnonparfun, x0);
        
        x_scale     = x_char/char_base;
        x_mean      = mean_base*x_scale;
        
        X.mean      = x_mean;
        X.scale     = x_scale;
        X.shift     = 0;
    case 33
        mean_base   = X.mean_base;
        char_base   = X.char_base;
        
        x_scale     = x_char/char_base;
        x_mean      = mean_base*x_scale;
        
        X.mean      = x_mean;
        X.scale     = x_scale;
        X.shift     = 0;
    otherwise
        error(['Unkown or not implemented distributiont type (dist): ', num2str(dist)])
end

% -------------------------------------------------------------------------
% NESTED functions
% -------------------------------------------------------------------------
    function P = my_lognormcdf(X, meanX, covX)
        mu_lognorm      = log(meanX.^2./sqrt((meanX.*covX).^2 + meanX.^2));
        sigma_lognorm   = sqrt(log((meanX.*covX).^2./meanX.^2 + 1));
        P               = logncdf(X, mu_lognorm, sigma_lognorm);
    end

    function P = my_gumbelcdf(x, meanX, covX)
        gamma           = 0.5772156649015328606065120900824024310421;
        beta            = sqrt(6)./pi*covX.*meanX;
        mupar           = meanX - beta.*gamma;
        P               = exp(-exp(-(x-mupar)./beta));
    end

    function f = my_parnonparfun(x)
       shift = 0;
       scale = x;
       f = nonparametric_cdf(x_char, dist_ID, shift, scale) - P_char;
    end
%     % maybe it could be reduced to fzero
%     function f = my_parnonparfun(x)
%        shift = x(1);
%        scale = x(2);
%        
%        % match the fractile with the non-exceedance probability
%        f(1) = nonparametric_cdf(x_frac, dist_ID, shift, scale) - P_frac;
%        % ensure a fixed cov
%        f(2) = scale*1*cov/(scale*1 + shift) - cov;   
%     end
end