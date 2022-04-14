% EC2 based shear resistance formula (proposal Yuguang 2019-October-11)
%
% VR = EC2_proposed_Yuguang_2019(fc, Asl, b, d, C_c, gamma_C)
%
% MIND THE UNITS! The formula is dimensionally inconsistent.
%
%INPUT
% fc        concrete compressive strength, [MPa]
% Asl       area of tensile reinforcement in the considered section, [mm^2]
% b         width, [mm]
% d         depth, [mm]
% C_c       model uncertainty factor, [-]
% gamma_C   partial factor, [-]
%
%OUTPUT
% VR    shear resistance, [kN]


function VR = EC2_proposed_Yuguang_2019(fc, Asl, b, d, C_c, gamma_C)

% % ad_ratio        = 3.0;  % treat as input par., later
D_lower         = 16;   % treat as input par., later

% normal weight concrete is assumed
d_dg            = 32*ones(size(fc));

idx             = fc > 60;
d_dg(idx)       = 16 + D_lower.*(60./fc(idx)).^(1/2);
d_dg            = bsxfun(@min, d_dg, 40);


% effective shear span
% % a_cs            = ad_ratio*d;
% % a_v             = reshape(max([a_cs(:), 2.5*d(:)],[],2), size(d));
   
rho_l           = Asl ./ (b .* d); 
% % rho_l(rho_l > 0.04) = 0.04; 

% calculate design shear resistance in [kN], according to eq. (6.4):
VR              = 1e-3 .* (C_c/gamma_C .* ((100 .* rho_l .* fc .* d_dg./d).^(1/3)) .* b .* d);


end