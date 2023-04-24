function [y, Ebank, Tthaw, Tent, Fthaw] = RunPfModel_NoSlump(H, U, Tw, S, zmax, dz, tmax, dt, f_ice, rhoB, Ti, Tau_crit, M)
%% run timeseries to calculate the rate of permafrost riverbank erosion by thaw- and entrainment-limited conditions
% inputs:
%   H = water depth (m)
%   U = average water velocity (m/s)
%   Tw = water temperature (degC)
%   S = channel slope (m/m)
%   zmax = bank height (m)
%   dz = vertical resolution for bank profile (m)
%   tmax = number of timesteps (days)
%   dt = timestep (days)
%   f_ice = mass fraction of riverbank comprised of ice (kg/kg)
%   rhoB = permafrost bank bulk density (kg/m3)
%   Ti = initial permafrost temperature (degC)
%   Tau_crit = critical shear stress to entrain sediment (Pa)
%   M = coefficient for entrainment equation (kg/m2/s)
% outputs:
%   y = matrix of updated bank ycoordinates (m) [zmax/dz x tmax/dt] 
%   Ebank = bank erosion rate (m/day)
%   Tthaw = fraction of time with thaw-limited erosion (days/days)
%   Tent = fraction of time with entrainment-limited erosion (days/days)
%   Fthaw = fraction of erosion that was thaw-limited (m/m)

f_sed = 1 - f_ice;                      % mass fraction of sediment in bank (kg/kg)
Lf = 334000 * f_ice;                    % permafrost latent heat of fusion (J/kg)
cp = f_ice*2108 + f_sed*700;            % specific heat of permafrost (J/kg/K)

y = zeros(zmax/dz+1, tmax/dt+1);
Ebank = zeros(zmax/dz+1, tmax/dt+1);
thawTF = zeros(1, tmax/dt);
for ts = 2:dt*tmax
    [y(:,ts), Ebank(:,ts), thawTF(ts-1)] = ...
        RunTimestep_NoSlump(y(:,ts-1),zmax,dz,Tw(ts),U(ts),H(ts),S,dt,f_sed,rhoB,Lf,cp,Ti,Tau_crit,M);
end

% calculate fraction of time thaw- vs entrainment-limited
Tthaw = sum(thawTF)/length(thawTF);
Tent = 1 - Tthaw;

% calculate bank fraction of erosion thaw- vs entrainment-limited
Fthaw = sum(thawTF.*max(Ebank(:,2:end)))/sum(max(Ebank(:,2:end)));

end