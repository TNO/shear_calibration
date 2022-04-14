% Eurocode advanced load combination rule
% EN1990:E:2002 6.4.3.2 Eq (6.10a-b)

function E = advanced_load_comb(gamma_G, ksi, G, gamma_Q1, psi01, Q1, gamma_Q2, psi02, Q2)

% -----------------------------
% Initilalize
% -----------------------------
if nargin < 9
    Q2 = 0;
end
if nargin < 8
    psi02 = 0;
end
if nargin < 7
    gamma_Q2 = 0;
end

% -----------------------------
% Calculate
% -----------------------------
% (6.10.a)
E1  = gamma_G.*G + gamma_Q1.*psi01.*Q1 + gamma_Q2.*psi02.*Q2;
%E1  = gamma_G.*G + gamma_Q1.*psi01.*Q1; + gamma_Q2.*psi02.*Q2;

% (6.10.b)
E21 = gamma_G.*ksi.*G + gamma_Q1.*Q1 + gamma_Q2.*psi02.*Q2; 
E22 = gamma_G.*ksi.*G + gamma_Q2.*Q2 + gamma_Q1.*psi01.*Q1;
%E21 = gamma_G.*ksi.*G + gamma_Q1.*Q1; + gamma_Q2.*psi02.*Q2; 
%E22 = gamma_G.*ksi.*G + gamma_Q2.*Q2; + gamma_Q1.*psi01.*Q1;
E2  = max(E21, E22);

E   = max(E1, E2);

end