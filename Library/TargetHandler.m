function [XtPPsMCI, COEtPPs, MEEtPPs, COEtdotsPPs] = TargetHandler(Xt_MCI, COEt, MEEt, tspan, ...
          EarthPPsMCI, SunPPsMCI, MoonPPsECI, deltaE, psiM, deltaM, muE, muS)


% Recall Global Variables
global DU TU muM Rm

% Initialize Local Variables
COEtdots = zeros(length(tspan), 6);
omegas_LVLH = zeros(length(tspan), 3);

% Compute Local Variables
for i = 1 : length(tspan)

    ti = tspan(i);      % get current epoch

    % Compute Third, Fourth Body and Moon Harmonics Perturbing Accelerations
    a34B = ThirdFourthBody(MEEt(i, :), ti, EarthPPsMCI, SunPPsMCI, muE, muS);
    aG_M = MoonHarmPerts(MEEt(i, :), MoonPPsECI, ti, muM, deltaE, psiM, deltaM);
    ap = a34B + aG_M;

    % Retrieve COE Values
    [~, ~, incl, ~, omega, nu] = S2C(COEt);
    theta_t = omega + nu;

    % Compute COE Derivatives
    COEtdots = get_COEdots(COEt(i, :), Xt_MCI(i, :), ap);
    [a_dot, e_dot, incl_dot, Omega_dot, omega_dot, nu_dot] = S2C(COEtdots);
    theta_t_dot = omega_dot + nu_dot;

    % Compute LVLH Angular Velocity wrt MCI
    omega_LVLH = [Omega_dot*sin(incl)*sin(theta_t) + incl_dot*cos(theta_tt); ...
                  Omega_dot*sin(incl)*cos(theta_t) - incl_dot*sin(theta_t); ...
                  Omega_dot*cos(incl) + theta_t_dot];

    % Store Necessary Quantities
    COEtdots(i, :) = [a_dot, e_dot, incl_dot, Omega_dot, omega_dot, nu_dot];
    omegas_LVLH(i, :) = omega_LVLH';

end

% Perform the Interpolation
XtPPsMCI = get_statePP(tspan, Xt_MCI);
COEtPPs = get_statePP(tspan, COEt);
MEEtPPs = get_statePP(tspan, MEEt);
COEtdotsPPs = get_statePP(tspan, COEtdots);
omegaPPsLVLH = get_statePP(tspan, omegas_LVLH);

% Compute the Derivative of LVLH's Angular Velocity wrt MCI
% !!!


end