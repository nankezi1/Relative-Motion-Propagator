function [R, Rdot] = get_R_Rdot(Xt_MCI, t, EarthPPsMCI, SunPPsMCI, MoonPPsECI, muE, muS, deltaE, psiM, deltaM)
% Description: this function retrieves the Rotation Matrix from MCI to LVLH
% and its derivative Rdot.
% 
% Inputs:
% Xt_MCI = target carthesian state in MCI
% t = current epoch
% EarthPPsMCI = pp struct with the interpolation of the Earth State in MCI
% SunPPsMCI = pp struct with the interpolation of the Sun State in MCI
% MoonPPsECI = pp struct with the interpolation of the Moon State in ECI
% muE = Earth's gravitational parameter in canonical units
% muS = Sun's gravitational parameter in canonical units
% deltaE = Earth's Ecliptic Obliquity
% psiM = Moon's Precession Angle
% deltaM = Moon's Equator Obliquity
% 
% Outputs:
% R = MCI2LVLH rotation matrix
% R_dot = MCI2LVLH rotation matrix derivative

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

% fprintf('The norm of R is: %.3f\nThe norm of Rdot is: %.3f\n\n', norm(R), norm(Rdot))


R_dot = [sin(Omega)*sin(incl)*sin(theta_t)*incl_dot- cos(Omega)*sin(theta_t)*theta_t_dot - cos(Omega)*cos(incl)*sin(theta_t)*Omega_dot - sin(Omega)*cos(incl)*cos(theta_t)*theta_t_dot - sin(Omega)*cos(theta_t)*Omega_dot,   cos(Omega)*cos(theta_t)*Omega_dot - sin(Omega)*sin(theta_t)*theta_t_dot - sin(Omega)*cos(incl)*sin(theta_t)*Omega_dot + cos(Omega)*cos(incl)*cos(theta_t)*theta_t_dot - cos(Omega)*sin(incl)*sin(theta_t)*incl_dot, cos(incl)*sin(theta_t)*incl_dot + cos(theta_t)*sin(incl)*theta_t_dot
         sin(Omega)*sin(theta_t)*Omega_dot - cos(Omega)*cos(theta_t)*theta_t_dot - cos(Omega)*cos(incl)*cos(theta_t)*Omega_dot + sin(Omega)*cos(theta_t)*sin(incl)*incl_dot + sin(Omega)*cos(incl)*sin(theta_t)*theta_t_dot, - cos(Omega)*sin(theta_t)*Omega_dot - sin(Omega)*cos(theta_t)*theta_t_dot - sin(Omega)*cos(incl)*cos(theta_t)*Omega_dot - cos(Omega)*cos(theta_t)*sin(incl)*incl_dot - cos(Omega)*cos(incl)*sin(theta_t)*theta_t_dot, cos(incl)*cos(theta_t)*incl_dot - sin(incl)*sin(theta_t)*theta_t_dot
                                                                                                                                                                                      cos(Omega)*sin(incl)*Omega_dot + sin(Omega)*cos(incl)*incl_dot,                                                                                                                                                                                         sin(Omega)*sin(incl)*Omega_dot - cos(Omega)*cos(incl)*incl_dot, -sin(incl)*incl_dot];



end