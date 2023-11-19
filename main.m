%% Relative Motion Propagator - Leonardo Russo

close all
clear
clc

addpath('Library/')
addpath('Data/')
addpath('Data/Planets/')

skip = 0;

if skip == 0

%% Hyperparameters and Settings

% Define options for ode113()
OptionsODE = odeset('RelTol', 1e-7, 'AbsTol', 1e-6, 'MaxStep', Inf);

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
Npoints = 10000;
Nperiods = 1.2;         % to set the final time
Nperiods = 0.5;

% Interpolate the Ephemeris and Retrieve Target's Initial State
[X0t_MCI, COE0t, MEE0t, EarthPPsMCI, DSGPPsMCI, SunPPsMCI, MoonPPsECI, time, t0, tf] = ...
 EphemerisHandler(deltaE, psiM, deltaM, Npoints, Nperiods);


% Introduce the Chaser Trajectory
% RHO0_MCI = random_delta(-10, 10, -1e-3, 1e-3);          % km, km/s
RHO0_MCI = [1 0 0 0 0 0]';                              % km, km/s
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
[XtPPsMCI, COEtPPs, MEEtPPs, COEtdotsPPs, omegaPPsLVLH] = ...
    TargetHandler(Xt_MCI, COEt, MEEt, tspan, EarthPPsMCI, SunPPsMCI, MoonPPsECI, deltaE, psiM, deltaM, muE, muS);

% Convert initial conditions into LVLH
RHO0_MCI = X0c_MCI - X0t_MCI;
RHO0_LVLH = rhoMCI2LVLH(RHO0_MCI, X0t_MCI, t0, EarthPPsMCI, SunPPsMCI, MoonPPsECI, muE, muS, deltaE, psiM, deltaM);


% Create Progress Bar for Chaser Propagation
global pbar
pbar = waitbar(0, 'Performing the Chaser Trajectory Propagation');

% Perform the Integration
[~, RHO_LVLH] = ode113(@(t, RHO_LVLH) DynamicalModelRM(t, RHO_LVLH, EarthPPsMCI, SunPPsMCI, ...
                                      muE, muS, tspan, MoonPPsECI, deltaE, ...
                                      psiM, deltaM, t0, tf, XtPPsMCI, omegaPPsLVLH), tspan, RHO0_LVLH, OptionsODE);
close(pbar);


% Conversion of the Results
Xc_MCI = zeros(size(RHO_LVLH, 1), 6);
RHO_MCI = zeros(size(RHO_LVLH, 1), 6);
omega_LVLH = zeros(size(RHO_LVLH, 1), 3);
omegadot_LVLH = zeros(size(RHO_LVLH, 1), 3);

for i = 1 : size(RHO_LVLH, 1)
    
    RHO_MCI(i, :) = rhoLVLH2MCI(RHO_LVLH(i, :)', Xt_MCI(i, :)', tspan(i), EarthPPsMCI, SunPPsMCI, MoonPPsECI, muE, muS, deltaE, psiM, deltaM)';
    Xc_MCI(i, :) = Xt_MCI(i, :) + RHO_MCI(i, :);

    % Compute Angular Velocity of LVLH wrt MCI and its derivative
    [omega_r, omegadot_r] = PolyEval(tspan(i), tspan, flip(omegaPPsLVLH(1).coefs, 2));
    [omega_t, omegadot_t] = PolyEval(tspan(i), tspan, flip(omegaPPsLVLH(2).coefs, 2));
    [omega_h, omegadot_h] = PolyEval(tspan(i), tspan, flip(omegaPPsLVLH(3).coefs, 2));
    omega_LVLH(i, :) = [omega_r, omega_t, omega_h]';
    omegadot_LVLH(i, :) = [omegadot_r, omegadot_t, omegadot_h]';

end

RHOref_MCI = Xrefc_MCI - Xt_MCI;


%% Visualization of the Results

close all
clc

% Draw the Target, Chaser and Reference Chaser Trajectories in MCI
figure('name', 'Trajectory in MCI Space')
T = DrawTrajMCI3D(Xt_MCI(:, 1:3)*DU, '#d1d1d1', '-.');
Cref = DrawTrajMCI3D(Xrefc_MCI(:, 1:3)*DU, '#61f4ff');
C = DrawTrajMCI3D(Xc_MCI(:, 1:3)*DU);
legend([T, Cref, C], {'Target Trajectory', 'Chaser Reference Trajectory', 'Chaser Trajectory'}, 'location', 'best');
if savechoice
    saveas(gcf, strcat('Output/TrajectoriesMCI.jpg'))
end


% Plot the evolution of the Chaser State in LVLH
figure('name', 'Chaser Trajectory in LVLH Space')
C_LVLH = DrawTrajLVLH3D(RHO_LVLH(:, 1:3)*DU);
if savechoice
    saveas(gcf, strcat('Output/TrajectoryLVLH.jpg'))
end


% Visualize LVLH State
figure('name', 'Chaser LVLH State')
subplot(2, 1, 1)
plot((tspan - t0) * TU / Day, RHO_LVLH(:, 1)*DU)
hold on
grid on
plot((tspan - t0) * TU / Day, RHO_LVLH(:, 2)*DU)
plot((tspan - t0) * TU / Day, RHO_LVLH(:, 3)*DU)
xlabel('t [Days]')
ylabel('[km]')
legend('$\delta_r$', '$\delta_\theta$', '$\delta_h$', 'interpreter', 'latex', 'fontsize', 12, 'location', 'best')

subplot(2, 1, 2)
plot((tspan - t0) * TU / Day, RHO_LVLH(:, 4)*DU/TU)
hold on
grid on
plot((tspan - t0) * TU / Day, RHO_LVLH(:, 5)*DU/TU)
plot((tspan - t0) * TU / Day, RHO_LVLH(:, 6)*DU/TU)
xlabel('t [Days]')
ylabel('[km/s]')
legend('$\delta_{\dot{r}}$', '$\delta_{\dot{\theta}}$', '$\delta_{\dot{h}}$', 'interpreter', 'latex', 'fontsize', 12, 'location', 'best')
if savechoice
    saveas(gcf, strcat('Output/RHO State LVLH.jpg'))
end


% Visualize MCI State Error
figure('name', 'Chaser MCI State Error wrt MEE Reference Propagation')
subplot(2, 1, 1)
plot((tspan - t0) * TU / Day, (Xc_MCI(:, 1) - Xrefc_MCI(:, 1)) * DU)
hold on
grid on
plot((tspan - t0) * TU / Day, (Xc_MCI(:, 2) - Xrefc_MCI(:, 2)) * DU)
plot((tspan - t0) * TU / Day, (Xc_MCI(:, 3) - Xrefc_MCI(:, 3)) * DU)
xlabel('t [Days]')
ylabel('[km]')
legend('$\delta_x$', '$\delta_y$', '$\delta_z$', 'interpreter', 'latex', 'fontsize', 12, 'location', 'best')

subplot(2, 1, 2)
plot((tspan - t0) * TU / Day, (Xc_MCI(:, 4) - Xrefc_MCI(:, 4)) * DU/TU)
hold on
grid on
plot((tspan - t0) * TU / Day, (Xc_MCI(:, 5) - Xrefc_MCI(:, 5)) * DU/TU)
plot((tspan - t0) * TU / Day, (Xc_MCI(:, 6) - Xrefc_MCI(:, 6)) * DU/TU)
xlabel('t [Days]')
ylabel('[km/s]')
legend('$\delta_{v_x}$', '$\delta_{v_y}$', '$\delta_{v_z}$', 'interpreter', 'latex', 'fontsize', 12, 'location', 'best')
if savechoice
    saveas(gcf, strcat('Output/State MCI Error.jpg'))
end


% Visualize RHO_MCI wrt MCI State Error
figure('name', 'RHO_MCI State wrt RHOref_MCI State')
subplot(2, 1, 1)
plot((tspan - t0) * TU / Day, (RHO_MCI(:, 1) - RHOref_MCI(:, 1))*DU)
hold on
grid on
plot((tspan - t0) * TU / Day, (RHO_MCI(:, 2) - RHOref_MCI(:, 2))*DU)
plot((tspan - t0) * TU / Day, (RHO_MCI(:, 3) - RHOref_MCI(:, 3))*DU)
xlabel('t [Days]')
ylabel('[km]')
legend('$\delta_x$', '$\delta_y$', '$\delta_z$', 'interpreter', 'latex', 'fontsize', 12, 'location', 'best')

subplot(2, 1, 2)
plot((tspan - t0) * TU / Day, (RHO_MCI(:, 4) - RHOref_MCI(:, 4))*DU/TU)
hold on
grid on
plot((tspan - t0) * TU / Day, (RHO_MCI(:, 5) - RHOref_MCI(:, 5))*DU/TU)
plot((tspan - t0) * TU / Day, (RHO_MCI(:, 6) - RHOref_MCI(:, 6))*DU/TU)
xlabel('t [Days]')
ylabel('[km/s]')
legend('$\delta_{v_x}$', '$\delta_{v_y}$', '$\delta_{v_z}$', 'interpreter', 'latex', 'fontsize', 12, 'location', 'best')
if savechoice
    saveas(gcf, strcat('Output/RHO State MCI Error.jpg'))
end



% % Show the Evolution of omega_LVLH
% figure('name', 'Evolution of omega_LVLH')
% plot((tspan - tspan(1))*TU/Day, omega_LVLH/TU)
% xlabel('$t \ [days]$', 'Interpreter','latex', 'FontSize', 12)
% ylabel('$\omega^{(LVLH)} \ [rad/s]$', 'Interpreter','latex', 'FontSize', 12)
% legend('$\omega_{r}$', '$\omega_{\theta}$', '$\omega_{h}$', 'location', 'best', 'interpreter', 'latex', 'fontsize', 12)


% % Show the Evolution of omegadot_LVLH
% figure('name', 'Evolution of omegadot_LVLH')
% plot((tspan - tspan(1))*TU/Day, omegadot_LVLH/TU^2)
% xlabel('$t \ [days]$', 'Interpreter','latex', 'FontSize', 12)
% ylabel('$\dot{\omega}^{(LVLH)} \ [rad/s^2]$', 'Interpreter','latex', 'FontSize', 12)
% legend('$\dot{\omega}_{r}$', '$\dot{\omega}_{\theta}$', '$\dot{\omega}_{h}$', 'location', 'best', 'interpreter', 'latex', 'fontsize', 12)
