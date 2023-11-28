%% Relative Motion Propagator - Leonardo Russo

close all
clear
clc

addpath('Library/')
addpath('Data/')
addpath('Data/Planets/')
addpath('Data/Materials/')

% % To Do's
% - fix the dynamical model


%% Hyperparameters and Settings

% Define options for ode113()
OptionsODE = odeset('RelTol', 1e-7, 'AbsTol', 1e-6, 'MaxStep', Inf);

savechoice = 0;     % set as 1 to save a copy of the plots locally

% Define Global Variables
global DU TU Rm muM pbar log

% Define HyperParameters
DU = 10000;                                             % km
TU = sqrt(1 / 4902.7779 * DU^3);                        % s
Rm = 1738.1 / DU;                                       % DU
muM = 4902.7779 * TU^2 / DU^3;                          % DU^3/TU^2
muE = 398600.4418 * TU^2 / DU^3;                        % DU^3/TU^2
muS = 132712440018 * TU^2 / DU^3;                       % DU^3/TU^2
deltaE = deg2rad(23.45);                                % rad
psiM = deg2rad(-81.7 + 360/18.6 * (5 + 4.5/12));        % rad
deltaM = deg2rad(1.5);                                  % rad

Day = 86400;                                            % s


% Define the n° of points for the Interpolation
Npoints = 1000;
Nperiods = 1.2;         % to set the final time

% Interpolate the Ephemeris and Retrieve Target's Initial State
[X0t_MCI, COE0t, MEE0t, EarthPPsMCI, DSGPPsMCI, SunPPsMCI, MoonPPsECI, time, t0, tf, Npoints] = ...
 EphemerisHandler(deltaE, psiM, deltaM, Npoints, Nperiods);


% Define the Chaser Initial Conditions
RHO0_LVLH = [5e-2, 0, 0, 0, 0, 0]';      % km, km/s
RHO0_LVLH = [RHO0_LVLH(1:3)/DU; RHO0_LVLH(4:6)/DU*TU];

RHO0_MCI = rhoLVLH2MCI(RHO0_LVLH, X0t_MCI, t0, EarthPPsMCI, SunPPsMCI, MoonPPsECI, muE, muS, deltaE, psiM, deltaM);

X0c_MCI = X0t_MCI + RHO0_MCI;

COE0c = rvPCI2COE(X0c_MCI', muM)';

MEE0c = COE2MEE(COE0c')';


%% Propagate the Target and Chaser Trajectory using MEE

% Define the timespan for the propagation
tspan = linspace(t0, tf, Npoints);

% Create Progress Bar for Target Propagation
pbar = waitbar(0, 'Performing the Target Trajectory Propagation');

% Perform the Integration
[~, MEEt] = ode113(@(t, MEE) DynamicalModelMEE(t, MEE, EarthPPsMCI, SunPPsMCI, ...
                                               muE, muS, MoonPPsECI, deltaE, psiM, deltaM, ...
                                               t0, tf), tspan, MEE0t, OptionsODE);
close(pbar);

% Create Progress Bar for Chaser Propagation
pbar = waitbar(0, 'Performing the Chaser Reference Trajectory Propagation');

% Perform the Integration
[~, MEErefc] = ode113(@(t, MEE) DynamicalModelMEE(t, MEE, EarthPPsMCI, SunPPsMCI, ...
                                               muE, muS, MoonPPsECI, deltaE, psiM, deltaM, ...
                                               t0, tf), tspan, MEE0c, OptionsODE);
close(pbar);

% Conversion from MEE to COE
COEt = MEE2COE(MEEt);
COErefc = MEE2COE(MEErefc);

% Conversion from COE to MCI
Xt_MCI = COE2rvPCI(COEt, muM);
Xrefc_MCI = COE2rvPCI(COErefc, muM);

RHOref_MCI = Xrefc_MCI - Xt_MCI;


% % Show the Target and Chaser Trajectories in MCI
% figure('name', 'Trajectory in MCI Space')
% T = DrawTrajMCI3D(Xt_MCI(:, 1:3)*DU, '#d1d1d1', '-.');
% Cref = DrawTrajMCI3D(Xrefc_MCI(:, 1:3)*DU, '#61f4ff');
% legend([T, Cref], {'Target Trajectory', 'Chaser Reference Trajectory'}, 'location', 'best');
% view([140, 30]);

% save('Data/WorkspaceMEE.mat');


%% Propagate Chaser Trajectory using Relative Motion

% close all
% clear
% clc
% 
% load('Data/WorkspaceMEE.mat');

% Interpolate the Target Trajectory
[omegaPPsLVLH, omegadotPPsLVLH] = TargetHandler(Xt_MCI, COEt, MEEt, tspan, EarthPPsMCI, ...
                                                SunPPsMCI, MoonPPsECI, deltaE, psiM, deltaM, muE, muS);

% Create Progress Bar for Chaser Propagation
pbar = waitbar(0, 'Performing the Chaser Trajectory Propagation');
log = fopen('log.txt', 'w+');

% Combine the Target and Chaser States into Y
Y0 = [MEE0t; RHO0_LVLH];

% Perform the Integration
[~, Y] = ode113(@(t, Y) DynamicalModelRM(t, Y, EarthPPsMCI, SunPPsMCI, ...
                                         muE, muS, MoonPPsECI, deltaE, psiM, deltaM, ...
                                         t0, tf, omegaPPsLVLH, omegadotPPsLVLH), tspan, Y0, OptionsODE);

close(pbar);
fclose(log);

% Retrieve the Target and Chaser States
MEEt = Y(:, 1:6);
COEt = MEE2COE(MEEt);
Xt_MCI1 = COE2rvPCI(COEt, muM);
RHO_LVLH = Y(:, 7:12);

% Conversion of the Results
Xc_MCI = zeros(size(RHO_LVLH, 1), 6);
RHO_MCI = zeros(size(RHO_LVLH, 1), 6);
omega_LVLH = zeros(size(RHO_LVLH, 1), 3);
omegadot_LVLH = zeros(size(RHO_LVLH, 1), 3);

pbar = waitbar(0, 'Performing the Final Post-Processing Steps');
for i = 1 : size(RHO_LVLH, 1)
    
    RHO_MCI(i, :) = rhoLVLH2MCI(RHO_LVLH(i, :)', Xt_MCI(i, :)', tspan(i), EarthPPsMCI, SunPPsMCI, MoonPPsECI, muE, muS, deltaE, psiM, deltaM)';
    Xc_MCI(i, :) = Xt_MCI(i, :) + RHO_MCI(i, :);

    % Compute Angular Velocity of LVLH wrt MCI and its derivative
    omega_LVLH(i, :) = ppsval(omegaPPsLVLH, tspan(i))';
    omegadot_LVLH(i, :) = ppsval(omegadotPPsLVLH, tspan(i))';

    waitbarMessage = sprintf('Final Post-Processing Progress: %.2f%%\n', i/length(tspan)*100);
    waitbar(i/length(tspan), pbar, waitbarMessage);      % update the waitbar

end
close(pbar)


save('Data/WorkspaceRM.mat');


%% Visualization of the Results

load('Data/WorkspaceRM.mat');

close all

% Draw the Target, Chaser and Reference Chaser Trajectories in MCI
figure('name', 'Trajectory in MCI Space')
T = DrawTrajMCI3D(Xt_MCI(:, 1:3)*DU, '#d1d1d1', '-.');
Cref = DrawTrajMCI3D(Xrefc_MCI(:, 1:3)*DU, '#61f4ff');
C = DrawTrajMCI3D(Xc_MCI(:, 1:3)*DU);
legend([T, Cref, C], {'Target Trajectory', 'Chaser Reference Trajectory', 'Chaser Trajectory'}, 'location', 'best');
view([140, 30]);
if savechoice
    saveas(gcf, strcat('Output/Trajectories MCI.jpg'))
end


% Plot the evolution of the Chaser State in LVLH
figure('name', 'Chaser Trajectory in LVLH Space')
C_LVLH = DrawTrajLVLH3D(RHO_LVLH(:, 1:3)*DU);
if savechoice
    saveas(gcf, strcat('Output/Trajectory LVLH.jpg'))
end


% Visualize LVLH State
figure('name', 'Chaser LVLH State')
subplot(2, 1, 1)
plot((tspan - t0) * TU / Day, RHO_LVLH(:, 1)*DU)
hold on
grid on
plot((tspan - t0) * TU / Day, RHO_LVLH(:, 2)*DU)
plot((tspan - t0) * TU / Day, RHO_LVLH(:, 3)*DU)
title('Chaser LVLH Position')
xlabel('$t \ [days]$', 'interpreter', 'latex', 'fontsize', 12)
ylabel('$[km]$', 'interpreter', 'latex', 'fontsize', 12)
legend('$r$', '$\theta$', '$h$', 'interpreter', 'latex', 'fontsize', 12, 'location', 'best')

subplot(2, 1, 2)
plot((tspan - t0) * TU / Day, RHO_LVLH(:, 4)*DU/TU)
hold on
grid on
plot((tspan - t0) * TU / Day, RHO_LVLH(:, 5)*DU/TU)
plot((tspan - t0) * TU / Day, RHO_LVLH(:, 6)*DU/TU)
title('Chaser LVLH Velocity')
xlabel('$t \ [days]$', 'interpreter', 'latex', 'fontsize', 12)
ylabel('$[km/s]$', 'interpreter', 'latex', 'fontsize', 12)
legend('$\dot{r}$', '$\dot{\theta}$', '$\dot{h}$', 'interpreter', 'latex', 'fontsize', 12, 'location', 'best')
if savechoice
    saveas(gcf, strcat('Output/RHO State LVLH.jpg'))
end


% Visualize Chaser MCI State Error - Relative Motion vs MEE Propagation
figure('name', 'Chaser MCI State Error - Relative Motion vs MEE Propagation')
subplot(2, 1, 1)
plot((tspan - t0) * TU / Day, (Xc_MCI(:, 1) - Xrefc_MCI(:, 1)) * DU)
hold on
grid on
plot((tspan - t0) * TU / Day, (Xc_MCI(:, 2) - Xrefc_MCI(:, 2)) * DU)
plot((tspan - t0) * TU / Day, (Xc_MCI(:, 3) - Xrefc_MCI(:, 3)) * DU)
title('Chaser MCI Position Error')
xlabel('$t \ [days]$', 'interpreter', 'latex', 'fontsize', 12)
ylabel('$[km]$', 'interpreter', 'latex', 'fontsize', 12)
legend('$\delta_x$', '$\delta_y$', '$\delta_z$', 'interpreter', 'latex', 'fontsize', 12, 'location', 'best')

subplot(2, 1, 2)
plot((tspan - t0) * TU / Day, (Xc_MCI(:, 4) - Xrefc_MCI(:, 4)) * DU/TU)
hold on
grid on
plot((tspan - t0) * TU / Day, (Xc_MCI(:, 5) - Xrefc_MCI(:, 5)) * DU/TU)
plot((tspan - t0) * TU / Day, (Xc_MCI(:, 6) - Xrefc_MCI(:, 6)) * DU/TU)
title('Chaser MCI Velocity Error')
xlabel('$t \ [days]$', 'interpreter', 'latex', 'fontsize', 12)
ylabel('$[km/s]$', 'interpreter', 'latex', 'fontsize', 12)
legend('$\delta_{v_x}$', '$\delta_{v_y}$', '$\delta_{v_z}$', 'interpreter', 'latex', 'fontsize', 12, 'location', 'best')
if savechoice
    saveas(gcf, strcat('Output/Chaser MCI State Error - RM vs MEE Propagations.jpg'))
end


% Visualize the Evolution of omega_LVLH
figure('name', 'Evolution of omega_LVLH')
plot((tspan - tspan(1))*TU/Day, omega_LVLH/TU)
title('Evolution of omega_LVLH')
xlabel('$t \ [days]$', 'Interpreter','latex', 'FontSize', 12)
ylabel('$\omega^{(LVLH)} \ [rad/s]$', 'Interpreter','latex', 'FontSize', 12)
legend('$\omega_{r}$', '$\omega_{\theta}$', '$\omega_{h}$', 'location', 'best', 'interpreter', 'latex', 'fontsize', 12)
if savechoice
    saveas(gcf, strcat('Output/omega.jpg'))
end


% Visualize the Evolution of omegadot_LVLH
figure('name', 'Evolution of omegadot_LVLH')
plot((tspan - tspan(1))*TU/Day, omegadot_LVLH/TU^2)
title('Evolution of omegadot_LVLH')
xlabel('$t \ [days]$', 'Interpreter','latex', 'FontSize', 12)
ylabel('$\dot{\omega}^{(LVLH)} \ [rad/s^2]$', 'Interpreter','latex', 'FontSize', 12)
legend('$\dot{\omega}_{r}$', '$\dot{\omega}_{\theta}$', '$\dot{\omega}_{h}$', 'location', 'best', 'interpreter', 'latex', 'fontsize', 12)
if savechoice
    saveas(gcf, strcat('Output/omegadot.jpg'))
end


% Visualize RHO_MCI wrt reference RHO_MCI
figure('name', 'RHO_MCI State wrt RHOref_MCI State')
subplot(2, 1, 1)
plot((tspan - t0) * TU / Day, (RHO_MCI(:, 1) - RHOref_MCI(:, 1))*DU)
hold on
grid on
plot((tspan - t0) * TU / Day, (RHO_MCI(:, 2) - RHOref_MCI(:, 2))*DU)
plot((tspan - t0) * TU / Day, (RHO_MCI(:, 3) - RHOref_MCI(:, 3))*DU)
xlabel('$t \ [days]$', 'interpreter', 'latex', 'fontsize', 12)
ylabel('$[km]$', 'interpreter', 'latex', 'fontsize', 12)
legend('$\delta_x$', '$\delta_y$', '$\delta_z$', 'interpreter', 'latex', 'fontsize', 12, 'location', 'best')

subplot(2, 1, 2)
plot((tspan - t0) * TU / Day, (RHO_MCI(:, 4) - RHOref_MCI(:, 4))*DU/TU)
hold on
grid on
plot((tspan - t0) * TU / Day, (RHO_MCI(:, 5) - RHOref_MCI(:, 5))*DU/TU)
plot((tspan - t0) * TU / Day, (RHO_MCI(:, 6) - RHOref_MCI(:, 6))*DU/TU)
xlabel('$t \ [days]$', 'interpreter', 'latex', 'fontsize', 12)
ylabel('$[km]$', 'interpreter', 'latex', 'fontsize', 12)
legend('$\delta_{v_x}$', '$\delta_{v_y}$', '$\delta_{v_z}$', 'interpreter', 'latex', 'fontsize', 12, 'location', 'best')
if savechoice
    saveas(gcf, strcat('Output/RHO_MCI vs RHOref_MCI Error.jpg'))
end


%% Testing

load('Data/WorkspaceRM.mat');

testing = 1;

if testing

% Dynamical Model Tests
clc
M = size(Y, 1);
dMEEt = zeros(M, 6);
dRHO_LVLH = zeros(M, 6);

pbar = waitbar(0, 'Testing the Propagation');
for i = 1 : M

    t = tspan(i);
    
    % Retrieve Data from Input
    MEEt_i = MEEt(i, :)';
    RHO_LVLH_i = RHO_LVLH(i, :)';
    
    % Clock for the Integration
    Day = 86400;  % seconds in a day
    Hour = 3600;  % seconds in an hour
    Min = 60;     % seconds in a minute
    tDAY = floor((t - t0) * TU / Day);      % calculate the elapsed time components
    tHR = floor(((t - t0) * TU - tDAY * Day) / Hour);
    tMIN = floor(((t - t0) * TU - tDAY * Day - tHR * Hour) / Min);
    timeStr = sprintf('Time Elapsed: %02d days, %02d hrs, %02d mins', tDAY, tHR, tMIN);     % create a string for the time
    waitbarMessage = sprintf('Testing Progress: %.2f%%\n%s', (t-t0)/(tf-t0)*100, timeStr);      % create the waitbar message including the time and progress percentage
    waitbar((t-t0)/(tf-t0), pbar, waitbarMessage);      % update the waitbar
    
    % Retrieve RHO State Variables
    rho_LVLH_i = RHO_LVLH_i(1:3);
    rhodot_LVLH_i = RHO_LVLH_i(4:6);
    
    % Retrieve Target State in MCI
    COEt_i = MEE2COE(MEEt_i')';
    Xt_MCI_i = COE2rvPCI(COEt_i', muM)';
    rt_MCI_i = Xt_MCI_i(1:3);
    
    % Convert RHO state into MCI
    RHO_MCI_i = rhoLVLH2MCI(RHO_LVLH_i, Xt_MCI_i, t, EarthPPsMCI, SunPPsMCI, MoonPPsECI, muE, muS, deltaE, psiM, deltaM);
    
    % Compute Chaser State in MCI
    Xc_MCI_i = Xt_MCI_i + RHO_MCI_i;
    COEc_i = rvPCI2COE(Xc_MCI_i', muM)';
    MEEc_i = COE2MEE(COEc_i')';
    rc_MCI_i = Xc_MCI_i(1:3);
    
    % Target's Third, Fourth Body and Moon Harmonics Perturbing Accelerations
    a34Bt_i = ThirdFourthBody(MEEt_i, t, EarthPPsMCI, SunPPsMCI, muE, muS);
    aG_Mt_i = MoonHarmPerts(MEEt_i, MoonPPsECI, t, muM, deltaE, psiM, deltaM);
    apt_LVLH_i = a34Bt_i + aG_Mt_i;
    
    % Chaser's Third, Fourth Body and Moon Harmonics Perturbing Accelerations
    a34Bc_i = ThirdFourthBody(MEEc_i, t, EarthPPsMCI, SunPPsMCI, muE, muS);
    aG_Mc_i = MoonHarmPerts(MEEc_i, MoonPPsECI, t, muM, deltaE, psiM, deltaM);
    apc_LVLHc_i = a34Bc_i + aG_Mc_i;
    
    % Convert Perturbating Accelerations into MCI
    [R_MCI2LVLHt_i, ~] = get_R_Rdot(Xt_MCI_i, t, EarthPPsMCI, SunPPsMCI, MoonPPsECI, muE, muS, deltaE, psiM, deltaM);
    [R_MCI2LVLHc_i, ~] = get_R_Rdot(Xc_MCI_i, t, EarthPPsMCI, SunPPsMCI, MoonPPsECI, muE, muS, deltaE, psiM, deltaM);
    
    apc_MCI_i = R_MCI2LVLHc_i'*apc_LVLHc_i;
    
    % Compute Angular Velocity of LVLH wrt MCI and its derivative
    omega_LVLH_i = ppsval(omegaPPsLVLH, t);
    omegadot_LVLH_i = ppsval(omegadotPPsLVLH, t);
    
    % Rotate the Necessary Vectors
    rt_LVLH_i = R_MCI2LVLHt_i * rt_MCI_i;     % this is -r_t2M^(LVLH)
    apc_LVLHt_i = R_MCI2LVLHt_i * apc_MCI_i;
    
    % MEE Propagation Quantities
    eta_i = get_eta(MEEt_i);
    G_i = get_G(MEEt_i, muM, eta_i);
    
    % RM Propagation Quantities
    rt_i = norm(rt_LVLH_i);
    rc_i = norm(rc_MCI_i);
    q_i = (dot(rho_LVLH_i, rho_LVLH_i) + 2*dot(rho_LVLH_i, rt_LVLH_i))/rt_i^2;
    
    
    % Assign State Derivatives
    dMEEt(i, 1:6) = G_i*apt_LVLH_i;
    dMEEt(i, 6) = dMEEt(6) + sqrt(muM/MEEt_i(1)^3)*eta_i^2;
    dRHO_LVLH(i, 1:3) = rhodot_LVLH_i;
    dRHO_LVLH(i, 4:6) = -2*cross(omega_LVLH_i, rhodot_LVLH_i) - cross(omegadot_LVLH_i, rho_LVLH_i) - cross(omega_LVLH_i, cross(omega_LVLH_i, rho_LVLH_i)) + muM/rt_i^3*((q_i*(2+q_i+(1+q_i)^(1/2)))/((1+q_i)^(3/2)*((1+q_i)^(1/2)+1)))*rt_LVLH_i - muM/rc_i^3*rho_LVLH_i + apc_LVLHt_i - apt_LVLH_i;

end
close(pbar)

figure('name', 'Evolution of adimensionalized RHO State Derivatives')
plot((tspan - t0) * TU / Day, dRHO_LVLH(:, 1:3))
hold on
plot((tspan - t0) * TU / Day, dRHO_LVLH(:, 4:6))
title('Evolution of adimensionalized RHO State Derivatives')
xlabel('$t \ [days]$', 'interpreter', 'latex', 'fontsize', 12)
% ylabel('$[km/s]\, , \,[km/s^2]$', 'interpreter', 'latex', 'fontsize', 12)
legend('$\dot{r}$', '$\dot{\theta}$', '$\dot{h}$', '$\ddot{r}$', '$\ddot{\theta}$', '$\ddot{h}$', 'interpreter', 'latex', 'fontsize', 12, 'location', 'best')

if savechoice
    saveas(gcf, strcat('Output/Relative Motion Adimensionalized Derivatives.jpg'))
end

end

% %% Manual Testing
% 
% dRHO_LVLH_PPs = get_statePP(tspan, dRHO_LVLH);
% ddRHO_LVLHr_PPs = fnder(dRHO_LVLH_PPs(1), 1);
% ddRHO_LVLHt_PPs = fnder(dRHO_LVLH_PPs(2), 1);
% ddRHO_LVLHh_PPs = fnder(dRHO_LVLH_PPs(3), 1);
% ddRHO_LVLH_PPs = [ddRHO_LVLHr_PPs; ddRHO_LVLHt_PPs; ddRHO_LVLHh_PPs];
% 
% ddRHO_LVLH = zeros(length(tspan), 3);
% for i = 1 : length(tspan)
%     ddRHO_LVLH(i, :) = ppsval(ddRHO_LVLH_PPs, tspan(i))';
% end
% 
% figure('name', 'Evolution of RHO State Derivatives from Interpolation')
% plot((tspan - t0) * TU / Day, dRHO_LVLH(:, 1:3)*DU/TU)
% hold on
% plot((tspan - t0) * TU / Day, ddRHO_LVLH(:, 1:3)*DU/TU^2)
% title('Evolution of RHO State Derivatives from Interpolation')
% xlabel('$t \ [days]$', 'interpreter', 'latex', 'fontsize', 12)
% ylabel('$[km/s]\, , \,[km/s^2]$', 'interpreter', 'latex', 'fontsize', 12)
% legend('$\dot{r}$', '$\dot{\theta}$', '$\dot{h}$', '$\ddot{r}$', '$\ddot{\theta}$', '$\ddot{h}$', 'interpreter', 'latex', 'fontsize', 12, 'location', 'best')
