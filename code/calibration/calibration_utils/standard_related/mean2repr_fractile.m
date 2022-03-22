% Calculate the representative value from the mean and probabilistic and
% standardized models where the representative value is a fractile.

function [x_repr, X] = mean2repr_fractile(x_mean, X)

dist = X.dist;

switch dist
    case {1, 2, 3, 11}
        [x_repr, X] = mean2repr_fractile_parametric(x_mean, X);
    case {33}
        [x_repr, X] = mean2repr_fractile_nonparametric(x_mean, X);
    otherwise
        error(['Unkown or not implemented distributiont type (dist): ', num2str(dist)])
end

end