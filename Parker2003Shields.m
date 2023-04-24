function Shields_crit = Parker2003Shields(D50)
%% calculate the critical Shields stress to entrain non-cohesive sediment using the Parker et al. (2003) version of Brownlie (1981)
% input:
%   D50 = sediment median grain size (m)
% output:
%   Shields_crit = critical Shields stress (dimensionless)

% constants
g = 9.81;       % gravitational acceleration (m/s2)
R = 1.65;       % sediment submerged specific gravity (dimensionless)
nu = 1e-6;       % water kinematic viscosity (m/s2)

% calculate particle Reynolds number
Re_p = sqrt(R*g*D50^3)/nu;
% calculate critical Shields stress
Shields_crit = 0.5*(0.22*Re_p^(-0.6) + 0.06*10^(-7.7*Re_p^(-0.6)));

end