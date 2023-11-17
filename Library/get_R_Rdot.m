function [R, Rdot] = get_R_Rdot(Xt_MCI, t, EarthPPsMCI, SunPPsMCI, MoonPPsECI, muE, muS, deltaE, psiM, deltaM)
% Description: this function retrieves the Rotation Matrix from MCI to LVLH
% and its derivative Rdot.

% Recall Global Variables
global muM Rm DU TU

% Retrieve Data from Input
COEt = rvPCI2COE(Xt_MCI', muM)';
MEEt = COE2MEE(COEt')';
[~, ~, incl, Omega, omega, nu] = S2C(COEt);
theta_t = omega + nu;

% Compute Third, Fourth Body and Moon Harmonics Perturbing Accelerations
a34B = ThirdFourthBody(MEEt, t, EarthPPsMCI, SunPPsMCI, muE, muS);
aG_M = MoonHarmPerts(MEEt, MoonPPsECI, t, muM, deltaE, psiM, deltaM);
ap = a34B + aG_M;

% Compute COE Derivatives
COEtdots = get_COEdots(COEt, Xt_MCI, ap);
[~, ~, incl_dot, Omega_dot, omega_dot, nu_dot] = S2C(COEtdots);
theta_t_dot = omega_dot + nu_dot;

% Compute MCI2LVLH Rotation Matrix and Rotation Matrix Derivative
R = R3(theta_t) * R1(incl) * R3(Omega);
Rdot = R3dot(theta_t)*R1(incl)*R3(Omega)*theta_t_dot + ...
       R3(theta_t)*R1dot(incl)*R3(Omega)*incl_dot + ...
       R3(theta_t)*R1(incl)*R3dot(Omega)*Omega_dot;

end