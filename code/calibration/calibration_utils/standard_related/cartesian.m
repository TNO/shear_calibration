% Cartesian product of `n` vectors.
%
% [C C_mx] = CARTESIAN(x1, x2, ..., xn)
%
%INPUT
% x1, x2, ..., xn   vectors (1d array)
%
%OUTPUT
% C     all n-element combinations of the inputs, the columns are follow 
%       the order of the input vectors (2D matrix).
% C_mx  has the same elements as `C` but shaped into a matrix of `size`:
%       `(length(x1), length(x2), ..., length(xn), n)`. Without loss of
%       generality, for `n=2` `C(3, 2, :)` is the vector that is formed by
%       the third element of `x1` and the second element of `x2`.        
%
% Based on: https://stackoverflow.com/a/33283459/4063376

function [C, C_mx] = cartesian(varargin)
    args = varargin;
    n = nargin;

    [F{1:n}] = ndgrid(args{:});
    C_mx = cat(n+1, F{:});
    C = reshape(C_mx, [], n);
end