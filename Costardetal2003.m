function E = Costardetal2003(H, U, Tw, rhoB, Lf, cp, Ti)
%% calculate thaw-limited permafrost riverbank erosion rates using the formulation from Costard et al., (2003)
% inputs:
%   H = water depth (m)
%   U = average water velocity (m/s)
%   Tw = water temperature (degC)
%   rhoB = permafrost bank bulk density (kg/m3)
%   Lf = permafrost latent heat of fusion (J/kg)
%   cp = permafrost heat capacity (J/kg/degC)
%   Ti = initial permafrost temperature (degC)
% outputs:
%   E = thaw-limited erosion rate (m/s)

% dimensionless constants from Lunardini et al., 1986
A = 0.0078;
alpha = 0.3333;
beta = 0.9270;

% water material properties and constants
kappa_w = 0.6;          % thermal conductivity (W/m/degC)
nu = 1e-6;              % kinematic viscosity (m2/s)
Tf = 0;                 % heat of fusion (degC)
Pr = 10;                % Prandtl number, 7 - 13 for range of water temp

Re = U*H/nu;                                    % Reynolds number
Ch = A*Pr^(alpha)*Re^(beta)*kappa_w/H;          % heat transfer coefficient (W/m2/degC)
E = Ch*(Tw - Tf) / (rhoB * (Lf + cp*(Ti - Tf)));        % thaw-limited erosion rate (m/s)
if E < 0
    E = 0;
end

end