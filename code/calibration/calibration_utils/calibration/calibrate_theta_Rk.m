% Calibrate the parameters in the semi-probabilistic format
%
%SYNOPSYS
% [parm_hat, objfun_val, exitflag] = CALIBRATE_theta_Rk(Prob, Options)
%

function [calibr_par, objfun_val, exitflag] = calibrate_theta_Rk(Prob, Prob_actions, DS, Options)

verbose = Options.verbose;


if verbose == 0
    opt = optimoptions('fmincon', 'Display', 'final');
elseif verbose > 0
    opt = optimoptions('fmincon', 'Display', 'iter');
end

% this strict only to make comparing the impact of changes and verifying
% that the changes did not introduce errors
opt.StepTolerance = 1e-5;

x0              = 1.0;
lb              = 0.5;
ub              = 3;
[x, fval, exitflag]  = fmincon(@(x) obj_fun_theta_Rk(x, Prob, Prob_actions, DS, Options), x0, [],[],[],[],lb,ub,[],opt);

calibr_par      = x;
objfun_val      = fval;

end
