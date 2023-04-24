function [ynew, ynew_block, Ebank, Eblock, Eent, Ethaw, thawTF, bankTF] = RunTimestep_Vslump(yold, yold_block, ...
    zmax, dz, Tw, U, H, S, dt, f_sed, rhoB, Lf, cp, Ti, Tau_crit, M, sigma_T, sigma_C, sigma_S)
%% calculate the rate of permafrost riverbank erosion by thaw- and entrainment-limited conditions and considering slump blocks
% inputs:
%   yold = vector of initial bank ycoordinates (m) [n x 1] 
%   yold_block = initial slump block area (m2)
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
%   sigma_T = bank tensile strength (Pa)
%   sigma_C = bank compressive strength (Pa)
%   sigma_S = bank shear strength (Pa)
% outputs:
%   ynew = vector of updated bank ycoordinates (m) [n x 1] 
%   ynew_block = updated slump block area (m2)
%   Ebank = bank erosion rate (m/day)
%   Eblock = slump block erosion rate (m/day)
%   Eent = entrainment-limited erosion rate (m/day)
%   Ethaw = thaw-limited erosion rate (m/day)
%   thawTF = true/false (1/0) indicator if erosion was thaw-limited
%   bankTF = true/false (1/0) indicator if the bank was eroded, versus a slump block

% calculate erosion rates
time_corr = 60*60*24;                                           % time correction vector (s/day)
Ethaw = Costardetal2003(H,U,Tw,rhoB,Lf,cp,Ti)*time_corr;        % thaw-limited erosion rate (m/day)
Eent = Partheneides1965(H, S, Tau_crit, rhoB, f_sed, M)*time_corr;      % entrainment-limited erosion rate (m/day)

% determine if erosion thaw- or entrainment-limited
if Ethaw > Eent
    E = Eent*ones(size(yold));
    thawTF = 0;
else
    E = Ethaw*ones(size(yold));
    thawTF = 1;
end

% erode any slump blocks present on bank
if max(yold_block) > 0
    Eblock = mean(E)*dt;
    Ebank = 0;
    bankTF = 0;
    ynew_block = yold_block - H*Eblock;
    if min(ynew_block) < 0
        yextra = -ynew_block;
        ynew_block = 0;
        Ebank = yextra/H*ones(size(yold));
        Ebank(round(H/dz):end) = 0;
    end
    ynew = yold + Ebank;

% if no slump blocks, erode bank directly
else
    ynew_block = yold_block;
    Eblock = 0;
    if round(H/dz) < (zmax+dz)/dz
        E(round(H/dz):end) = 0;
    end
    Ebank = E*dt;
    bankTF = 1;
    ynew = yold + Ebank;
end

% Calculate the factor of safety for bank failure
[failTF, zfail, failAr] = BlockFailure(ynew, H, dz, zmax, rhoB, sigma_T, sigma_C, sigma_S);
% If bank failure occurs, remove that volume of sediment from the bank and
% store it as a slump block volume
if failTF == 1
    ynew_block = failAr;
    ynew(zfail+1:end) = ynew(zfail);
end

end