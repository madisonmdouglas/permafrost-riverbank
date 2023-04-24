function E = Partheneides1965(H, S, Tau_crit, rhoB, f_sed, M)
%% calculate entrainment-limited riverbank erosion for cohesive sediment under normal flow conditions using the entrainment equation from Partheneides (1965)
% inputs:
%   H = water depth (m)
%   S = channel slope (m/m)
%   Tau_crit = critical shear stress to entrain sediment (Pa)
%   rhoB = bulk density of riverbank (kg/m3)
%   f_sed = mass fraction of riverbank comprised of sediment (kg/kg)
%   M = coefficient for entrainment equation (kg/m2/s)
% outputs:
%   E = entrainment-limited bank erosion rate (m/s)

% constants
g = 9.81;           % gravitational acceleration (m/s2)
rho_w = 1000;       % water density (kg/m3)
n = 1;              % exponent for entrainment equation (dimensionless)

% calculate shear stress on the riverbank
epsilon = 0.2;
Tau_bank = rho_w*g*H*S/(1+epsilon);

% if bank stress is greater than the threshold for entrainment, calculate
% erosion rates
if Tau_bank > Tau_crit
    E = M*(Tau_bank / Tau_crit - 1)^n;
else
    E = 0;
end

% convert erosion rate from kg/m2/s to m/s
E = E / (rhoB*f_sed);

end