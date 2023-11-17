function [XtPPsMCI, COEtPPs, MEEtPPs, COEtdotsPPs] = TargetHandler(Xt_MCI, COEt, MEEt, tspan, ...
          EarthPPsMCI, SunPPsMCI, MoonPPsECI, deltaE, psiM, deltaM, muE, muS)


% Recall Global Variables
global DU TU muM Rm

% Initialize Local Variables
COEtdots = zeros(length(tspan), 6);

% Compute Local Variables
for i = 1 : length(tspan)

    ti = tspan(i);

    a = COEt(i, 1);
    e = COEt(i, 2);
    incl = COEt(i, 3);
    Omega = COEt(i, 4);
    omega = COEt(i, 5);
    nu = COEt(i, 6);
    theta_t = COEt(i, 5) + COEt(i, 6);
    
    rt_MCI = Xt_MCI(i, 1:3);
    vt_MCI = Xt_MCI(i, 4:6);
 
    pt = a*(1-e^2);
    rt = norm(rt_MCI);
    ht_vect = cross(rt_MCI, vt_MCI);
    ht = norm(ht_vect);

    % Compute Third, Fourth Body and Moon Harmonics Perturbing Accelerations
    a34B = ThirdFourthBody(MEEt(i, :), ti, EarthPPsMCI, SunPPsMCI, muE, muS);
    aG_M = MoonHarmPerts(MEEt(i, :), MoonPPsECI, ti, muM, deltaE, psiM, deltaM);
    
    ap = a34B + aG_M;
    ap_r = ap(1);
    ap_t = ap(2);
    ap_h = ap(3);

    % Compute COE Derivatives
    a_dot = 2*a^2/ht * (e*ap_r*sin(nu) + ap_t*(e*cos(nu) + 1));
    e_dot = sqrt(pt/muM)*ap_r*sin(nu) + sqrt(pt/muM)*ap_t*(e+e*cos(nu)^2+2*cos(nu))/(1+e*cos(nu));
    incl_dot = rt*ap_h/ht * cos(theta_t);
    Omega_dot = rt*ap_h/ht *sin(theta_t)/sin(incl);
    omega_dot = -rt*ap_h*sin(theta_t)*cos(incl)/(ht*sin(incl)) - ap_r/e*cos(nu)*sqrt(pt/muM) + ap_t*sqrt(pt/muM)*sin(nu)*(e*cos(nu)+2)/(e*(1+e*cos(nu)));
    nu_dot = sqrt(muM/pt^3)*(1+e*cos(nu))^2 + ap_r/e*cos(nu)*sqrt(pt/muM) - ap_t*sqrt(pt/muM)*sin(nu)*(e*cos(nu)+2)/e*(1+e*cos(nu));

    COEtdots(i, :) = [a_dot, e_dot, incl_dot, Omega_dot, omega_dot, nu_dot];

end

% Perform the Interpolation
XtPPsMCI = get_statePP(tspan, Xt_MCI);
COEtPPs = get_statePP(tspan, COEt);
MEEtPPs = get_statePP(tspan, MEEt);
COEtdotsPPs = get_statePP(tspan, COEtdots);

end