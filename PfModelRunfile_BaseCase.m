% Runfile for permafrost model hydraulics base case for sandy banks along
% the Yukon River at Stevens Village

% produces Figure 3

clear

% setup experimental run
rhoB = 861;                 % bulk density of mineral sediment (kg/m3), from Lininger et al. (2019)
rhoS = 2650;                % sediment density (kg/m3)
rhoW = 1000;                % water density (kg/m3)
g = 9.81;                   % gravitational acceleration (m/s2)
f_ice = 0.2362;             % mass fraction of ice for mineral sediment (kg/kg), from Lininger et al. (2019)
Ti = -1;                    % initial bank temperature (degC), from borehole data
S = 10.583e-5;              % channel slope (m/m), from Clement 1999 thesis

% sand
sigma_S = 5e4;                % riverbank shear strength (Pa)
sigma_C = 11.2e6;             % compressive strength (Pa), 1std = 4.1e6
sigma_T = 2.4e6;              % tensile strength (Pa), 1std = 0.2e6
D = 1e-3;                     % grain size (m)
Tau_crit = Parker2003Shields(D)*(rhoS - rhoW)*g*D;    % critical shear stress for entrainment (Pa)
M = 2.5e-5;                   % coefficient for entrainment equation (kg/m2/s)

% load hydraulic data for Yukon River at Stevens Village
load('StevensVillageforModel_v3.mat')
Tw = TempC;                         % water temperature (degC)
H = Depthm;                         % water depth (m)
Qw = DischargeMediancms;            % median water discharge (m3/s)
Qw25 = Discharge25cms;              % 25th percentile water discharge (m3/s)
Qw75 = Discharge75cms;              % 75th percentil water discharge (m3/s)
zmax = round(max(H),1);             % bankfull height (m)
dz = 0.001;                         % bank height grid (m)
U = Velocityms;                     % average water speed (m/s)
dt = mean(diff(Day));               % simulation timestep (days)
tmax = max(Day);                    % number of days to run simulation
clear Depthm DischargeMediancms Discharge25cms Discharge75cms TempC Velocityms

% run model for given conditions
[y, yblock, Ebank, Eblock, Eent, Ethaw, Tbank, Tblock, Tthaw, Tent, Fbank, Fthaw] = ...
    RunPfModel_Vslump(H, U, Tw, S, zmax, dz, tmax, dt, f_ice, rhoB, Ti, Tau_crit, M, sigma_T, sigma_C, sigma_S);
Ebank_tot = sum(y(:,end)*dz)/zmax;

% run model without entrainment-limited erosion
[y_NE,yblock_NE,~,~,~,~,~] = RunPfModel_NoEnt(H, U, Tw, zmax, dz, tmax, dt, f_ice, rhoB, Ti, sigma_T, sigma_C, sigma_S);
Ebank_NEtot = sum(y_NE(:,end-1)*dz)/zmax;

% run model without slump blocks
[y_NS,~,~,~,~] = RunPfModel_NoSlump(H, U, Tw, S, zmax, dz, tmax, dt, f_ice, rhoB, Ti, Tau_crit, M);
Ebank_NStot = sum(y_NS(:,end-1)*dz)/zmax;

% run model without entrainment-limited erosion or slump blocks
[y_NENS,~] = RunPfModel_NoEntNoSlump(H, U, Tw, zmax, dz, tmax, dt, f_ice, rhoB, Ti);
Ebank_NENStot = sum(y_NENS(:,end-1)*dz)/zmax;

% plot results in multi-panel figure
figure('Renderer', 'painters', 'Position', [10 10 1000 1000])

subplot(4,2,1)
shadedErrorBar(Day',Qw',[Qw75'-Qw'; Qw' - Qw25'])
ylabel('{\it Q_{w,SV}} (m^3/s)')
set(gca,'FontSize',16)
text(1,10.5,'a','FontSize',24,'Units','characters')
xlim([0,365])
box on

subplot(4,2,3)
plot(Tw,'k-','LineWidth',2)
ylabel('{\it T_{w,SV}} (\circC)')
set(gca,'FontSize',16)
text(1,10.5,'b','FontSize',24,'Units','characters')
xlim([0,365])
box on

Eent2 = Eent;
Eent2(Eblock>0) = NaN;
Eent2(Eent2<1e-2) = 1e-2;
Ethaw2 = Ethaw;
Ethaw2(Eblock>0) = NaN;
Ethaw2(Ethaw2<1e-2) = 1e-2;
Ebank2 = max(Ebank);
Ebank2(Eblock>0) = NaN;
Ebank2(Ebank2<1e-2) = 1e-2;
subplot(4,2,5)
hold on
plot(Eent2,'LineWidth',6,'Color',[243,124,30]/255)
plot(Ethaw2,'LineWidth',6,'Color',[30,144,208]/255)
plot(Ebank2,'k-','LineWidth',2)
% legend({'{\itE_{ent}}','{\itE_{thaw}}','{\itE_{bank}}'},'location','southwest')
ylabel('{\it E_{bank}} (m/day)')
set(gca,'FontSize',16,'YScale','log')
text(1,10.5,'c','FontSize',24,'Units','characters')
yticks([1e-2,1e-1,1,10,100,1e3])
yticklabels({'0','10^{-1}','10^0','10^1','10^2','10^3'})
ylim([1e-2,1e3])
xlim([0,365])
box on
hold off

subplot(4,2,7)
hold on
plot(Eent.*(Eblock>0),'LineWidth',6,'Color',[243,124,30]/255)
plot(Ethaw.*(Eblock>0),'LineWidth',6,'Color',[30,144,208]/255)
plot(Eblock,'k-','LineWidth',2)
legend({'{\itE_{ent}}','{\itE_{thaw}}','minimum {\itE}'},'location','southwest')
xlabel('Day of the year')
ylabel('{\it E_{block}} (m/day)')
set(gca,'FontSize',16,'YScale','log')
text(1,10.5,'d','FontSize',24,'Units','characters')
yticks([1e-2,1e-1,1,10,100,1e3])
yticklabels({'0','10^{-1}','10^0','10^1','10^2','10^3'})
ylim([1e-2,1e3])
xlim([0,365])
box on
hold off

doi = [100, 150, 200, 225, 250, 300];
newcolors = [0.95 0.95 0.95; 0.9 0.9 0.9; 0.85 0.85 0.85; 0.8 0.8 0.8; ...
    0.75 0.75 0.75; 0.7 0.7 0.7; 0.65 0.65 0.65];
subplot(4,2,2:2:8);
hold on
for i = 1:length(doi)
    fill([y(:,doi(i));12;12],[0:dz:zmax,zmax,0],newcolors(i,:))
    plot([y(round(H(doi(i))/dz),doi(i))-1,y(round(H(doi(i))/dz),doi(i))], ...
        [H(doi(i)),H(doi(i))],'LineWidth',2,'Color',[0,114,178]/255)
end
ylim([0,11])
xlim([-2,12]);
xlabel('Distance eroded (m)')
ylabel('Height above channel bed (m)')
set(gca,'FontSize',16)
text(1,60,'e','FontSize',24,'Units','characters')
text(3,10,'Top of the eroding bank','FontSize',12)
text(-1.25,5,'River channel','FontSize',12,'rotation',90)
text(0.5,5,'100 days','FontSize',12,'HorizontalAlignment','left','rotation',90)
text(2.7,1,'150 days','FontSize',12,'HorizontalAlignment','left','rotation',90)
text(6.8,1,'200 days','FontSize',12,'HorizontalAlignment','left','rotation',90)
text(8.5,1,'225 days','FontSize',12,'HorizontalAlignment','left','rotation',90)
text(10.1,1,'250 days','FontSize',12,'HorizontalAlignment','left','rotation',90)
text(11.3,1,'300 days','FontSize',12,'HorizontalAlignment','left','rotation',90)
box on
hold off

