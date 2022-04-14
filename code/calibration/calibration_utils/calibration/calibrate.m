% Calibrate the parameters in the semi-probabilistic format
%
%SYNOPSYS
% [parm_hat, objfun_val, exitflag] = CALIBRATE(Prob, Options)
%

function [calibr_par, objfun_val, exitflag, output] = calibrate(Prob, Prob_actions, DS, Options)

verbose = Options.verbose;


if verbose == 0
%     opt     = optimset('Display', 'notify');
    opt = optimoptions('fmincon', 'Display', 'final');
elseif verbose > 0
%     opt     = optimset('Display', 'iter');
    opt = optimoptions('fmincon', 'Display', 'iter');
end

% this strict only to make comparing the impact of changes and verifying
% that the changes did not introduce errors
opt.StepTolerance = 1e-5;

x0              = 1.3;
% we need a constrained optmization solver that is gradient free (Matlab has no
% such built-in function)
% x       = fminsearch(@(x) obj_fun(x, Prob, Options), x0, opt);
lb              = 0.5;
ub              = 3;
[x, fval, exitflag, output]  = fmincon(@(x) obj_fun(x, Prob, Prob_actions, DS, Options), x0, [],[],[],[],lb,ub,[],opt);

calibr_par      = x;
objfun_val      = fval;

end
