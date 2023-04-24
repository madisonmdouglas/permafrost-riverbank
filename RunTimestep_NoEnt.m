function [ynew, ynew_block, Ebank, Eblock, bankTF] = RunTimestep_NoEnt(yold, yold_block, ...
    zmax, dz, Tw, U, H, dt, rhoB, Lf, cp, Ti, sigma_T, sigma_C, sigma_S)
%% calculate the rate of permafrost riverbank erosion by thaw-limited conditions and considering slump blocks
% inputs:
%   yold = vector of initial bank ycoordinates (m) [n x 1] 
%   yold_block = initial slump block area (m2)
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
%   sigma_T = bank tensile strength (Pa)
%   sigma_C = bank compressive strength (Pa)
%   sigma_S = bank shear strength (Pa)
% outputs:
%   ynew = vector of updated bank ycoordinates (m) [n x 1] 
%   ynew_block = updated slump block area (m2)
%   Ebank = bank erosion rate (m/day)
%   Eblock = slump block erosion rate (m/day)
%   bankTF = true/false (1/0) indicator if the bank was eroded, versus a slump block

% calculate erosion rates
time_corr = 60*60*24*ones(size(yold));                          % time correction vector (s/day)
E = Costardetal2003(H,U,Tw,rhoB,Lf,cp,Ti)*time_corr;            % thaw-limited erosion rate (m/day)

% erode any slump blocks present on bank
if max(yold_block) > 0
    Eblock = mean(E)*dt;
    Ebank = 0;
    bankTF = 0;
    ynew_block = yold_block - H*Eblock;
    if min(ynew_block) < 0
        yextra = -ynew_block;
        ynew_block = 0;
        Ebank = yextra/H;
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