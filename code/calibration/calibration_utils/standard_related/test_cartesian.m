% Test our custom Cartestian product function (`cartesian`).
%
% Some of the test requires the Deep Learning Toolbox (function `combvec`).
clc; close all; clearvars

x = [2, 1];
y = [4, 5, 6, 3];
z = [11, 10, 12];

% ..........................................
% 2D
% ..........................................
[cp, cp_mx] = cartesian(x, y);
cp_expected = combvec(x, y)';
assert(all(all(cp == cp_expected)))

for ii = 1:length(x)
    for jj = 1:length(y)
        assert(all(squeeze(cp_mx(ii,jj,:))' == [x(ii), y(jj)]))
    end
end

% ..........................................
% 3D
% ..........................................
[cp, cp_mx] = cartesian(x, y, z);
cp_expected = combvec(x, y, z)';
assert(all(all(cp == cp_expected)))

for ii = 1:length(x)
    for jj = 1:length(y)
        for kk = 1:length(z)
            assert(all(squeeze(cp_mx(ii,jj,kk,:))' == [x(ii), y(jj), z(kk)]))
        end
    end
end





