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
% - understand yesterday's changes

skip = 0;

%% Hyperparameters and Settings

if skip == 0

% Define options for ode113()
OptionsODE = odeset('RelTol', 1e-10, 'AbsTol', 1e-9, 'MaxStep', Inf);

savechoice = 1;     % set as 1 to save a copy of the plots locally

% Define Global Variables
global DU TU Rm muM

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
Npoints = 100000;
Nperiods = 1.2;         % to set the final time

% Interpolate the Ephemeris and Retrieve Target's Initial State
[X0t_MCI, COE0t, MEE0t, EarthPPsMCI, DSGPPsMCI, SunPPsMCI, MoonPPsECI, time, t0, tf, Npoints] = ...
 EphemerisHandler(deltaE, psiM, deltaM, Npoints, Nperiods);


% Define the Chaser Initial Conditions
RHO0_LVLH = [5e-2, 0, 0, 0, 0, 0]';      % km, km/s
RHO0_LVLH = [RHO0_LVLH(1:3)/DU; RHO0_LVLH(4:6)/DU*TU];

RHO0_MCI = rhoLVLH2MCI(RHO0_LVLH, X0t_MCI, t0, EarthPPsMCI, SunPPsMCI, MoonPPsECI, muE, muS, deltaE, psiM, deltaM);

X0c_MCI = X0t_MCI + [RHO0_MCI(1:3)/DU; RHO0_MCI(4:6)/DU*TU];

COE0c = rvPCI2COE(X0c_MCI', muM)';

MEE0c = COE2MEE(COE0c')';


%% Propagate the Target and Chaser Trajectory using MEE

% Define the timespan for the propagation
tspan = linspace(t0, tf, Npoints);

% Create Progress Bar for Target Propagation
global pbar
pbar = waitbar(0, 'Performing the Target Trajectory Propagation');

% Perform the Integration
[~, MEEt] = ode113(@(t, MEE) DynamicalModelMEE(t, MEE, EarthPPsMCI, SunPPsMCI, ...
                                             muE, muS, time, MoonPPsECI, deltaE, ...
                                             psiM, deltaM, t0, tf), tspan, MEE0t, OptionsODE);
close(pbar);

% Create Progress Bar for Chaser Propagation
global pbar
pbar = waitbar(0, 'Performing the Chaser Reference Trajectory Propagation');

% Perform the Integration
[~, MEEc] = ode113(@(t, MEE) DynamicalModelMEE(t, MEE, EarthPPsMCI, SunPPsMCI, ...
                                             muE, muS, time, MoonPPsECI, deltaE, ...
                                             psiM, deltaM, t0, tf), tspan, MEE0c, OptionsODE);
close(pbar);

% Conversion from MEE to COE
COEt = MEE2COE(MEEt);
COErefc = MEE2COE(MEEc);

% Conversion from COE to MCI
Xt_MCI = COE2rvPCI(COEt, muM);
Xrefc_MCI = COE2rvPCI(COErefc, muM);


% Show the Target and Chaser Trajectories in MCI
figure('name', 'Trajectory in MCI Space')
T = DrawTrajMCI3D(Xt_MCI(:, 1:3)*DU, '#d1d1d1', '-.');
Cref = DrawTrajMCI3D(Xrefc_MCI(:, 1:3)*DU, '#61f4ff');
legend([T, Cref], {'Target Trajectory', 'Chaser Reference Trajectory'}, 'location', 'best');

save('Data/WorkspaceMEE.mat');

end

%% Propagate Chaser Trajectory using Relative Motion

close all
clear
clc

load('Data/WorkspaceMEE.mat');

% Interpolate the Target Trajectory
[XtPPsMCI, COEtPPs, MEEtPPs, COEtdotsPPs, omegaPPsLVLH, omegadotPPsLVLH] = ...
    TargetHandler(Xt_MCI, COEt, MEEt, tspan, EarthPPsMCI, SunPPsMCI, MoonPPsECI, deltaE, psiM, deltaM, muE, muS);

% Create Progress Bar for Chaser Propagation
global pbar log
pbar = waitbar(0, 'Performing the Chaser Trajectory Propagation');
log = fopen('log.txt', 'w+');

% Combine the Target and Chaser States into Y
Y0 = [MEE0t; RHO0_LVLH];

% Perform the Integration
[~, Y] = ode113(@(t, Y) DynamicalModelRM(t, Y, EarthPPsMCI, SunPPsMCI, ...
                                         muE, muS, tspan, MoonPPsECI, deltaE, ...
                                         psiM, deltaM, t0, tf, omegaPPsLVLH, omegadotPPsLVLH), tspan, Y0, OptionsODE);
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


%% Visualization of the Results

close all

% Draw the Target, Chaser and Reference Chaser Trajectories in MCI
figure('name', 'Trajectory in MCI Space')
T = DrawTrajMCI3D(Xt_MCI(:, 1:3)*DU, '#d1d1d1', '-.');
Cref = DrawTrajMCI3D(Xrefc_MCI(:, 1:3)*DU, '#61f4ff');
C = DrawTrajMCI3D(Xc_MCI(:, 1:3)*DU);
legend([T, Cref, C], {'Target Trajectory', 'Chaser Reference Trajectory', 'Chaser Trajectory'}, 'location', 'best');
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
xlabel('$t \ [days]$', 'interpreter', 'latex', 'fontsize', 12)
ylabel('$[km]$', 'interpreter', 'latex', 'fontsize', 12)
legend('$r$', '$\theta$', '$h$', 'interpreter', 'latex', 'fontsize', 12, 'location', 'best')

subplot(2, 1, 2)
plot((tspan - t0) * TU / Day, RHO_LVLH(:, 4)*DU/TU)
hold on
grid on
plot((tspan - t0) * TU / Day, RHO_LVLH(:, 5)*DU/TU)
plot((tspan - t0) * TU / Day, RHO_LVLH(:, 6)*DU/TU)
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
xlabel('$t \ [days]$', 'interpreter', 'latex', 'fontsize', 12)
ylabel('$[km]$', 'interpreter', 'latex', 'fontsize', 12)
legend('$\delta_x$', '$\delta_y$', '$\delta_z$', 'interpreter', 'latex', 'fontsize', 12, 'location', 'best')

subplot(2, 1, 2)
plot((tspan - t0) * TU / Day, (Xc_MCI(:, 4) - Xrefc_MCI(:, 4)) * DU/TU)
hold on
grid on
plot((tspan - t0) * TU / Day, (Xc_MCI(:, 5) - Xrefc_MCI(:, 5)) * DU/TU)
plot((tspan - t0) * TU / Day, (Xc_MCI(:, 6) - Xrefc_MCI(:, 6)) * DU/TU)
xlabel('$t \ [days]$', 'interpreter', 'latex', 'fontsize', 12)
ylabel('$[km/s]$', 'interpreter', 'latex', 'fontsize', 12)
legend('$\delta_{v_x}$', '$\delta_{v_y}$', '$\delta_{v_z}$', 'interpreter', 'latex', 'fontsize', 12, 'location', 'best')
if savechoice
    saveas(gcf, strcat('Output/Chaser MCI State Error - RM vs MEE Propagations.jpg'))
end


% Visualize the Evolution of omega_LVLH
figure('name', 'Evolution of omega_LVLH')
plot((tspan - tspan(1))*TU/Day, omega_LVLH/TU)
xlabel('$t \ [days]$', 'Interpreter','latex', 'FontSize', 12)
ylabel('$\omega^{(LVLH)} \ [rad/s]$', 'Interpreter','latex', 'FontSize', 12)
legend('$\omega_{r}$', '$\omega_{\theta}$', '$\omega_{h}$', 'location', 'best', 'interpreter', 'latex', 'fontsize', 12)
if savechoice
    saveas(gcf, strcat('Output/omega.jpg'))
end


% Visualize the Evolution of omegadot_LVLH
figure('name', 'Evolution of omegadot_LVLH')
plot((tspan - tspan(1))*TU/Day, omegadot_LVLH/TU^2)
xlabel('$t \ [days]$', 'Interpreter','latex', 'FontSize', 12)
ylabel('$\dot{\omega}^{(LVLH)} \ [rad/s^2]$', 'Interpreter','latex', 'FontSize', 12)
legend('$\dot{\omega}_{r}$', '$\dot{\omega}_{\theta}$', '$\dot{\omega}_{h}$', 'location', 'best', 'interpreter', 'latex', 'fontsize', 12)
if savechoice
    saveas(gcf, strcat('Output/omegadot.jpg'))
end


% % Visualize RHO_MCI wrt MCI State Error
% figure('name', 'RHO_MCI State wrt RHOref_MCI State')
% subplot(2, 1, 1)
% plot((tspan - t0) * TU / Day, (RHO_MCI(:, 1) - RHOref_MCI(:, 1))*DU)
% hold on
% grid on
% plot((tspan - t0) * TU / Day, (RHO_MCI(:, 2) - RHOref_MCI(:, 2))*DU)
% plot((tspan - t0) * TU / Day, (RHO_MCI(:, 3) - RHOref_MCI(:, 3))*DU)
% xlabel('$t \ [days]$', 'interpreter', 'latex', 'fontsize', 12)
% ylabel('$[km]$', 'interpreter', 'latex', 'fontsize', 12)
% legend('$\delta_x$', '$\delta_y$', '$\delta_z$', 'interpreter', 'latex', 'fontsize', 12, 'location', 'best')
% 
% subplot(2, 1, 2)
% plot((tspan - t0) * TU / Day, (RHO_MCI(:, 4) - RHOref_MCI(:, 4))*DU/TU)
% hold on
% grid on
% plot((tspan - t0) * TU / Day, (RHO_MCI(:, 5) - RHOref_MCI(:, 5))*DU/TU)
% plot((tspan - t0) * TU / Day, (RHO_MCI(:, 6) - RHOref_MCI(:, 6))*DU/TU)
% xlabel('$t \ [days]$', 'interpreter', 'latex', 'fontsize', 12)
% ylabel('$[km]$', 'interpreter', 'latex', 'fontsize', 12)
% legend('$\delta_{v_x}$', '$\delta_{v_y}$', '$\delta_{v_z}$', 'interpreter', 'latex', 'fontsize', 12, 'location', 'best')
% if savechoice
%     saveas(gcf, strcat('Output/RHO State MCI Error.jpg'))
% end


%% Testing

testing = 0;

if testing

% Dynamical Model Tests
clc
M = size(Xt_MCI, 1);
dRHO_LVLH = zeros(M, 6);

pbar = waitbar(0, 'Testing the Propagation');
for i = 1 : M

    t = tspan(i);
    
    % Clock for the Integration
    Day = 86400;  % seconds in a day
    Hour = 3600;  % seconds in an hour
    Min = 60;     % seconds in a minute
    tDAY = floor((t - t0) * TU / Day);      % calculate the elapsed time components
    tHR = floor(((t - t0) * TU - tDAY * Day) / Hour);
    tMIN = floor(((t - t0) * TU - tDAY * Day - tHR * Hour) / Min);
    timeStr = sprintf('Time Elapsed: %02d days, %02d hrs, %02d mins', tDAY, tHR, tMIN);     % create a string for the time
    waitbarMessage = sprintf('Progress: %.2f%%\n%s', (t-t0)/(tf-t0)*100, timeStr);      % create the waitbar message including the time and progress percentage
    waitbar((t-t0)/(tf-t0), pbar, waitbarMessage);      % update the waitbar
    
    % Retrieve RHO State Variables
    rho_LVLH = RHO_LVLH(i, 1:3)';
    rhodot_LVLH = RHO_LVLH(i, 4:6)';
    
    % Retrieve Target State in MCI
    Xt_MCIi = ppsval(XtPPsMCI, t);
    COEti = rvPCI2COE(Xt_MCIi', muM)';
    MEEti = COE2MEE(COEti')';
    rt_MCI = Xt_MCIi(1:3);
    rt = norm(rt_MCI);
    
    % Convert RHO state into MCI
    RHO_MCIi = rhoLVLH2MCI(RHO_LVLH(i, :)', Xt_MCIi, t, EarthPPsMCI, SunPPsMCI, MoonPPsECI, muE, muS, deltaE, psiM, deltaM);
    rho_MCI = RHO_MCIi(1:3);
    
    % Compute Chaser State in MCI
    Xc_MCIi = Xt_MCIi + RHO_MCIi;
    COEci = rvPCI2COE(Xc_MCIi', muM)';
    MEEci = COE2MEE(COEci')';
    rc_MCI = Xc_MCIi(1:3);
    rc = norm(rc_MCI);
    
    q = (dot(rho_MCI, rho_MCI) + 2*dot(rho_MCI, rt_MCI))/rt^2;
    
    % Target's Third, Fourth Body and Moon Harmonics Perturbing Accelerations
    a34Bt = ThirdFourthBody(MEEti, t, EarthPPsMCI, SunPPsMCI, muE, muS);
    aG_Mt = MoonHarmPerts(MEEti, MoonPPsECI, t, muM, deltaE, psiM, deltaM);
    apt_LVLH = a34Bt + aG_Mt;
    
    % Chaser's Third, Fourth Body and Moon Harmonics Perturbing Accelerations
    a34Bc = ThirdFourthBody(MEEci, t, EarthPPsMCI, SunPPsMCI, muE, muS);
    aG_Mc = MoonHarmPerts(MEEci, MoonPPsECI, t, muM, deltaE, psiM, deltaM);
    apc_LVLH = a34Bc + aG_Mc;
    
    % Convert Perturbating Accelerations into MCI
    [R_MCI2LVLH, ~] = get_R_Rdot(Xt_MCIi, t, EarthPPsMCI, SunPPsMCI, MoonPPsECI, muE, muS, deltaE, psiM, deltaM);
    apt_MCI = R_MCI2LVLH'*apt_LVLH;
    apc_MCI = R_MCI2LVLH'*apc_LVLH;
    
    % Compute Angular Velocity of LVLH wrt MCI and its derivative
    [omega_r, omegadot_r] = PolyEval(t, tspan, flip(omegaPPsLVLH(1).coefs, 2));
    [omega_t, omegadot_t] = PolyEval(t, tspan, flip(omegaPPsLVLH(2).coefs, 2));
    [omega_h, omegadot_h] = PolyEval(t, tspan, flip(omegaPPsLVLH(3).coefs, 2));
    omega_LVLHi = [omega_r, omega_t, omega_h]';
    omegadot_LVLH = [omegadot_r, omegadot_t, omegadot_h]';
    
    % Assign State Derivatives
    dRHO_LVLH(i, 1:3) = rhodot_LVLH;
    dRHO_LVLH(i, 4:6) = -2*cross(omega_LVLHi, rhodot_LVLH) - cross(omegadot_LVLH, rho_LVLH) - cross(omega_LVLHi, cross(omega_LVLHi, rho_LVLH)) + ...
                muM/rt^3*((q*(2+q+(1+q)^(1/2)))/((1+q)^(3/2)*((1+q)^(1/2)+1)))*rt_MCI - muM/rc^3*rho_MCI + apc_MCI - apt_MCI;

end
close(pbar)

figure('name', 'Evolution of RHO State Derivatives')
plot((tspan - t0) * TU / Day, dRHO_LVLH(:, 1:3)*DU)
hold on
plot((tspan - t0) * TU / Day, dRHO_LVLH(:, 4:6)*DU/TU)
xlabel('$t \ [days]$', 'interpreter', 'latex', 'fontsize', 12)
ylabel('$[km]\, , \,[km/s]$', 'interpreter', 'latex', 'fontsize', 12)
legend('$\dot{r}$', '$\dot{\theta}$', '$\dot{h}$', '$\ddot{r}$', '$\ddot{\theta}$', '$\ddot{h}$', 'interpreter', 'latex', 'fontsize', 12, 'location', 'best')

end