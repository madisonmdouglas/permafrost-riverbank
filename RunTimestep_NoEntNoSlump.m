function [ynew, Ebank] = RunTimestep_NoEntNoSlump(yold, zmax, dz, Tw, U, H, dt, rhoB, Lf, cp, Ti)
%% calculate the rate of permafrost riverbank erosion by thaw-limited conditions
% inputs:
%   yold = vector of initial bank ycoordinates (m) [n x 1] 
%   zmax = bank height (m)
%   dz = vertical resolution for bank profile (m)
%   Tw = water temperature (degC)
%   U = average water velocity (m/s)
%   H = water depth (m)
%   dt = timestep (days)
%   rhoB = permafrost bank bulk density (kg/m3)
%   Lf = permafrost latent heat of fusion (J/kg)
%   cp = permafrost heat capacity (J/kg/degC)
%   Ti = initial permafrost temperature (degC)
% outputs:
%   ynew = vector of updated bank ycoordinates (m) [n x 1] 
%   Ebank = bank erosion rate (m/day)

% calculate erosion rates
time_corr = 60*60*24*ones(size(yold));                          % time correction vector (s/day)
E = Costardetal2003(H,U,Tw,rhoB,Lf,cp,Ti)*time_corr;            % thaw-limited erosion rate (m/day)

% erode riverbank
if round(H/dz) < (zmax+dz)/dz
    E(round(H/dz):end) = 0;
end
Ebank = E*dt;
ynew = yold + Ebank;

end