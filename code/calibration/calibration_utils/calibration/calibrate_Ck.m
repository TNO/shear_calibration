% Calibrate the parameters in the semi-probabilistic format
%
%SYNOPSYS
% [parm_hat, objfun_val, exitflag] = CALIBRATE_CK(Prob, Options)
%

function [calibr_par, objfun_val, exitflag] = calibrate_Ck(Prob, Prob_actions, DS, Options)

verbose = Options.verbose;


if verbose == 0
    opt = optimoptions('fmincon', 'Display', 'final');
elseif verbose > 0
    opt = optimoptions('fmincon', 'Display', 'iter');
end

x0              = 0.2;
lb              = 0.001;
ub              = 0.8;
[x, fval, exitflag]  = fmincon(@(x) obj_fun_Ck(x, Prob, Prob_actions, DS, Options), x0, [],[],[],[],lb,ub,[],opt);

calibr_par      = x;
objfun_val      = fval;

end

