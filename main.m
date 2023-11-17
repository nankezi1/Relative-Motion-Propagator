%% Relative Motion Propagator - Leonardo Russo

close all
clear
clc

addpath('Library/')
addpath('Data/')
addpath('Data/Planets/')

%% Hyperparameters and Settings

% Define options for ode45()
OptionsODE = odeset('RelTol', 1e-13, 'AbsTol', 1e-11, 'MaxStep', Inf);

savechoice = 0;     % set as 1 to save a copy of the plots locally

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
Nperiods = 2.2;         % to set the final time

% Interpolate the Ephemeris and Retrieve Target's Initial State
[X0t_MCI, COE0t, MEE0t, EarthPPsMCI, DSGPPsMCI, SunPPsMCI, MoonPPsECI, time, t0, tf] = ...
 EphemerisHandler(deltaE, psiM, deltaM, Npoints, Nperiods);


% Introduce the Chaser Trajectory
Rho0_MCI = random_delta(-10, 10, -1e-3, 1e-3);          % km, km/s
% rho0_MCI = [1 0 0 0 0 0]';        % km, km/s
X0c_MCI = X0t_MCI + [Rho0_MCI(1:3)/DU; Rho0_MCI(4:6)/DU*TU];

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
pbar = waitbar(0, 'Performing the Chaser Trajectory Propagation');

% Perform the Integration
[~, MEEc] = ode113(@(t, MEE) DynamicalModelMEE(t, MEE, EarthPPsMCI, SunPPsMCI, ...
                                             muE, muS, time, MoonPPsECI, deltaE, ...
                                             psiM, deltaM, t0, tf), tspan, MEE0c, OptionsODE);
close(pbar);

% Conversion from MEE to COE
COEt = MEE2COE(MEEt);
COEc = MEE2COE(MEEc);

% Conversion from COE to MCI
Xt_MCI = COE2rvPCI(COEt, muM);
Xc_MCI = COE2rvPCI(COEc, muM);


% Show the Target and Chaser Trajectories in MCI
figure('name', 'Trajectory in MCI Space')
T = DrawTrajMCI3D(Xt_MCI(:, 1:3)*DU, '#d1d1d1', '-.');
C_MEE = DrawTrajMCI3D(Xc_MCI(:, 1:3)*DU, '#61f4ff');
legend([T, C_MEE], {'Target Trajectory', 'Chaser Reference Trajectory'}, 'location', 'best');


%% Propagate Chaser Trajectory using Relative Motion

% Interpolate the Target Trajectory
TargetPPsMCI = get_statePP(tspan, Xt_MCI);

% Convert initial conditions into LVLH
Rho0_MCI = X0c_MCI - X0t_MCI;
Rho0_LVLH = rhoMCI2LVLH(Rho0_MCI, X0t_MCI, t0, EarthPPsMCI, SunPPsMCI, MoonPPsECI, muE, muS, deltaE, psiM, deltaM);


%% Perform the Integration

% [~, Rho_LVLH] = ode113(@(t, Rho_LVLH) DynamicalModelRM(t, Rho_LVLH, EarthPPsMCI, SunPPsMCI, ...
%                                       muE, muS, time, MoonPPsECI, deltaE, ...
%                                       psiM, deltaM, t0, tf, XtPPsMCI, COEtPPs, MEEtPPs, COEtdotsPPs), tspan, MEE0c, OptionsODE);
