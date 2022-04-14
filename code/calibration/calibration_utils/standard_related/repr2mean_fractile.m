% Calculate the mean from representative value and probabilistic and
% standardized models where the representative value is a fractile.

function [x_mean, X] = repr2mean_fractile(x_repr, X)

dist = X.dist;

switch dist
    case {1, 2, 3, 11}
        [x_mean, X] = repr2mean_fractile_parametric(x_repr, X);
    case {33}
        [x_mean, X] = repr2mean_fractile_nonparametric(x_repr, X);
    otherwise
        error(['Unkown or not implemented distributiont type (dist): ', num2str(dist)])
end

end