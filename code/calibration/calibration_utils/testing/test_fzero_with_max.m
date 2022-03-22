% test root finding

clearvars
close all
clc

x       = linspace(0,10,1e3);

y1      = 1 + 1*x;
y2      = -10 + 3*x;

% idx = y2 >= y1;
% y(idx) = y2(idx);

y       = max(y1, y2);

myfun   = @(xx) interp1(x, y, xx);

a1       = 3;
[xsol1, fsol1]      = fzero(@(xx) myfun(xx) - a1, 5);

a2       = 10;
[xsol2, fsol2]      = fzero(@(xx) myfun(xx) - a2, 5);

plot(x, y)
hold on
plot([min(x), max(x)], [a1, a1], 'red')
plot(xsol1, fsol1 + a1, 'ro')
hold on
plot([min(x), max(x)], [a2, a2], 'blue')
plot(xsol2, fsol2 + a2, 'bo')