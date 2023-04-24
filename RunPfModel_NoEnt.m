function [y, yblock, Ebank, Eblock, Tbank, Tblock, Fbank] = ...
    RunPfModel_NoEnt(H, U, Tw, zmax, dz, tmax, dt, f_ice, rhoB, Ti, sigma_T, sigma_C, sigma_S)
%% run timeseries to calculate the rate of permafrost riverbank erosion by thaw-limited conditions and considering slump blocks
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
%   sigma_T = bank tensile strength (Pa)
%   sigma_C = bank compressive strength (Pa)
%   sigma_S = bank shear strength (Pa)
% outputs:
%   y = matrix of updated bank ycoordinates (m) [zmax/dz x tmax/dt] 
%   y_block = vector of slump block area (m2) [1 x tmax/dt]
%   Ebank = bank erosion rate (m/day)
%   Eblock = slump block erosion rate (m/day)
%   Tbank = fraction of time spent eroding bank (days/days)
%   Tblock = fraction of time spent eroding slump blocks (days/days)
%   Fbank = fraction of erosion on the bank versus on blocks (m/m)

f_sed = 1 - f_ice;                      % mass fraction of sediment in bank (kg/kg)
Lf = 334000 * f_ice;                    % permafrost latent heat of fusion (J/kg)
cp = f_ice*2108 + f_sed*700;            % specific heat of permafrost (J/kg/K)

y = zeros(zmax/dz+1, tmax/dt+1);
yblock = zeros(1, tmax/dt+1);
Ebank = zeros(zmax/dz+1, tmax/dt+1);
Eblock = zeros(1, tmax/dt+1);
bankTF = zeros(1, tmax/dt);
for ts = 2:dt*tmax
    [y(:,ts), yblock(ts), Ebank(:,ts), Eblock(ts), bankTF(ts-1)] = ...
        RunTimestep_NoEnt(y(:,ts-1),yblock(ts-1),zmax,dz,Tw(ts),U(ts),H(ts),dt,rhoB,Lf,cp,Ti,sigma_T,sigma_C,sigma_S);
end

% calculate fraction of time spent eroding bank vs slump blocks
Tbank = sum(bankTF)/length(bankTF);
Tblock = 1 - Tbank;

% calculate bank fraction of erosion on the bank versus slump blocks
Fbank = sum(max(Ebank))/(sum(max(Ebank)) + sum(Eblock));

end