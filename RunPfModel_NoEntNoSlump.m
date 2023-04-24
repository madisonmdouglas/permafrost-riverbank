function [y, Ebank] = RunPfModel_NoEntNoSlump(H, U, Tw, zmax, dz, tmax, dt, f_ice, rhoB, Ti)
%% run timeseries to calculate the rate of permafrost riverbank erosion by thaw-limited conditions
% inputs:
%   H = water depth (m)
%   U = average water velocity (m/s)
%   Tw = water temperature (degC)
%   zmax = bank height (m)
%   dz = vertical resolution for bank profile (m)
%   tmax = number of timesteps (days)
%   dt = timestep (days)
%   f_ice = mass fraction of riverbank comprised of ice (kg/kg)
%   rhoB = permafrost bank bulk density (kg/m3)
%   Ti = initial permafrost temperature (degC)
% outputs:
%   y = matrix of updated bank ycoordinates (m) [zmax/dz x tmax/dt] 
%   Ebank = bank erosion rate (m/day)

f_sed = 1 - f_ice;                      % mass fraction of sediment in bank (kg/kg)
Lf = 334000 * f_ice;                    % permafrost latent heat of fusion (J/kg)
cp = f_ice*2108 + f_sed*700;            % specific heat of permafrost (J/kg/K)

y = zeros(zmax/dz+1, tmax/dt+1);
Ebank = zeros(zmax/dz+1, tmax/dt+1);
for ts = 2:dt*tmax
    [y(:,ts), Ebank(:,ts)] = RunTimestep_NoEntNoSlump(y(:,ts-1),zmax,dz,Tw(ts),U(ts),H(ts),dt,rhoB,Lf,cp,Ti);
end

end