function [ynew, Ebank, thawTF] = RunTimestep_NoSlump(yold, zmax, dz, Tw, U, H, S, dt, f_sed, rhoB, Lf, cp, Ti, Tau_crit, M)
%% calculate the rate of permafrost riverbank erosion by thaw- and entrainment-limited conditions
% inputs:
%   yold = vector of initial bank ycoordinates (m) [n x 1] 
%   zmax = bank height (m)
%   dz = vertical resolution for bank profile (m)
%   Tw = water temperature (degC)
%   U = average water velocity (m/s)
%   H = water depth (m)
%   S = channel slope (m/m)
%   dt = timestep (days)
%   f_sed = mass fraction of riverbank comprised of sediment (kg/kg)
%   rhoB = permafrost bank bulk density (kg/m3)
%   Lf = permafrost latent heat of fusion (J/kg)
%   cp = permafrost heat capacity (J/kg/degC)
%   Ti = initial permafrost temperature (degC)
%   Tau_crit = critical shear stress to entrain sediment (Pa)
%   M = coefficient for entrainment equation (kg/m2/s)
% outputs:
%   ynew = vector of updated bank ycoordinates (m) [n x 1] 
%   Ebank = bank erosion rate (m/day)
%   thawTF = true/false (1/0) indicator if erosion was thaw-limited

% calculate erosion rates
time_corr = 60*60*24*ones(size(yold));                          % time correction vector (s/day)
Ethaw = Costardetal2003(H,U,Tw,rhoB,Lf,cp,Ti)*time_corr;        % thaw-limited erosion rate (m/day)
Eent = Partheneides1965(H, S, Tau_crit, rhoB, f_sed, M)*time_corr;  % entrainment-limited erosion rate (m/day)

% determine if erosion thaw- or entrainment-limited
if Ethaw > Eent
    E = Eent;
    thawTF = 0;
else
    E = Ethaw;
    thawTF = 1;
end

% erode riverbank
if round(H/dz) < (zmax+dz)/dz
    E(round(H/dz):end) = 0;
end
Ebank = E*dt;
ynew = yold + Ebank;

end