% Eurocode simple load combination rule
% EN1990:E:2002 6.4.3.2 Eq (6.10)

function [E, E1, E2] = simple_load_comb(gamma_G, G, gamma_Q1, psi01, Q1, gamma_Q2, psi02, Q2)

% -----------------------------
% Initilalize
% -----------------------------
if nargin < 7
    Q2 = 0;
end
if nargin < 6
    psi02 = 0;
end
if nargin < 5
    gamma_Q2 = 0;
end

% -----------------------------
% Calculate
% -----------------------------
% (6.10)
E1 = gamma_G.*G + gamma_Q1.*Q1 + gamma_Q2.*psi02.*Q2;
E2 = gamma_G.*G + gamma_Q2.*Q2 + gamma_Q1.*psi01.*Q1;

E  = max(E1, E2);

end