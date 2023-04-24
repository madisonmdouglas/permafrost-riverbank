function [failTF, zfail, failAr] = BlockFailure(y, ztest, dz, zmax, rhoB, sigma_T, sigma_C, sigma_S)
%% calculate bank failure factor of safety and volume for irregular failure geometry based on Patsingahasanee et al. (2016)
% inputs:
%   yold = vector of bank ycoordinates (m) [n x 1] 
%   ztest = elevation of water surface which intersects failure plane (m)
%   dz = vertical resolution for bank profile (m)
%   zmax = bank height (m)
%   rhoB = permafrost bank bulk density (kg/m3)
%   sigma_T = bank tensile strength (Pa)
%   sigma_C = bank compressive strength (Pa)
%   sigma_S = bank shear strength (Pa)
% outputs:
%   failTF = true/false (1/0) if bank failure occurs
%   zfail = vector index for base of slump block
%   failAr = area of sediment that failed (m2)

% calculate the potential slump block geometry
zfail = round(ztest/dz);
if zfail < length(y)
    failBlock = y(zfail) - y(zfail+1:end);
    failAr = sum(failBlock)*dz;
else
    failAr = 0;
end

% screen for edge cases
if failAr <= 0
    failTF = 0;
    zfail = 0;
    failAr = 0;

% determine if a slump block occurs
else
    g = 9.81;           % gravitational acceleration (m/s2)

    % calculate center of mass coordinates
    zCOM = sum(dz*failBlock.*transpose(dz*zfail:dz:zmax)) / failAr;
    yCOM = sum(dz*failBlock.^2/2) / failAr;

    % calculate the distances to the point of rotation
    Lt = zmax - zCOM;
    Lc = zCOM - ztest;
    Lb = y(zfail) - yCOM;

    % calculate factor of safety for rotational failure
    Fs_rot = 2*(rhoB*g*failAr*Lb)/(sigma_C*Lc^2 + sigma_T*Lt^2);
    
    % calculate factor of safety for brittle failure
    Fs_brit =  rhoB*g*failAr / (sigma_S*(Lt+Lc));
    
    % if rotational or brittle failure occurs, block fails
    Fs = max([Fs_rot, Fs_brit]);
    failTF = (Fs>=1);
    
end

end

