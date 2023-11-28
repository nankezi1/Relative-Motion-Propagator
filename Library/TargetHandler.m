function [omegaPPsLVLH, omegadotPPsLVLH] = TargetHandler(Xt_MCI, COEt, MEEt, tspan, ...
          EarthPPsMCI, SunPPsMCI, MoonPPsECI, deltaE, psiM, deltaM, muE, muS)


% Recall Global Variables
global DU TU muM Rm

Day = 86400;    % s

% Initialize Local Variables
omega_LVLH = zeros(length(tspan), 3);

% Compute Local Variables
pbar = waitbar(0, 'Performing the Target Trajectory Interpolation');

for i = 1 : length(tspan)

    ti = tspan(i);      % get current epoch

    % Compute Third, Fourth Body and Moon Harmonics Perturbing Accelerations
    a34B = ThirdFourthBody(MEEt(i, :)', ti, EarthPPsMCI, SunPPsMCI, muE, muS);
    aG_M = MoonHarmPerts(MEEt(i, :)', MoonPPsECI, ti, muM, deltaE, psiM, deltaM);
    ap_LVLH = a34B + aG_M;

    % Retrieve COE Values
    [~, ~, incl, ~, omega, nu] = S2C(COEt);
    theta_t = omega + nu;

    % Compute COE Derivatives
    COEtdots = get_COEdots(COEt(i, :), Xt_MCI(i, :), ap_LVLH);
    [~, ~, incl_dot, Omega_dot, omega_dot, nu_dot] = S2C(COEtdots);
    theta_t_dot = omega_dot + nu_dot;

    % Compute LVLH Angular Velocity wrt MCI
    omega_LVLH(i, :) = [Omega_dot*sin(incl)*sin(theta_t) + incl_dot*cos(theta_t); ...
                        Omega_dot*sin(incl)*cos(theta_t) - incl_dot*sin(theta_t); ...
                        Omega_dot*cos(incl) + theta_t_dot]';

    waitbarMessage = sprintf('Target Interpolation Progress: %.2f%%\n', i/length(tspan)*100);
    waitbar(i/length(tspan), pbar, waitbarMessage);      % update the waitbar

end

% Perform the Interpolation
omegaPPsLVLH = get_statePP(tspan, omega_LVLH);

omegadot_rPPs = fnder(omegaPPsLVLH(1), 1);
omegadot_tPPs = fnder(omegaPPsLVLH(2), 1);
omegadot_hPPs = fnder(omegaPPsLVLH(3), 1);
omegadotPPsLVLH = [omegadot_rPPs; omegadot_tPPs; omegadot_hPPs];

close(pbar);

end