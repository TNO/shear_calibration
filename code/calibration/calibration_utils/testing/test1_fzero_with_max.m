% test root finding

clearvars
close all
clc

y1      = @(x) 1 + 1*x;
y2      = @(x) -3.5 + 3*x;

y       = @(x) max(y1(x), y2(x));

a1      = 0;
[xsol, fsol] = fzero(@(x) y(x) - a1, 5);
[xsol1, fsol1] = fzero(@(x) y1(x) - a1, 5);
[xsol2, fsol2] = fzero(@(x) y2(x) - a1, 5);


diff_y1 = y(xsol1) - a1
diff_y2 = y(xsol2) - a1


x       = linspace(0,1e1,1e3);
plot(x, y(x), 'k')
hold on
plot(x, y1(x), 'r:')
plot(x, y2(x), 'b:')
plot([xsol1, xsol1], [-10, 20], 'r')
plot([xsol2, xsol2], [-10, 20], 'b')
