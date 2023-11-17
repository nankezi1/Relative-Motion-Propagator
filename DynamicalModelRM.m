function dRho_LVLH = DynamicalModelRM(t, Rho_LVLH, EarthPPsMCI, SunPPsMCI, muE, ...
                                  muS, MoonPPsECI, deltaE, psiM, ...
                                  deltaM, t0, tf, XtPPsMCI, COEtPPs, MEEtPPs, COEtdotsPPs)


% Retrieve Global Variables
global muM Rm DU TU pbar

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


% Retrieve Rho State Variables
rho_LVLH = Rho_LVLH(1:3);
rhodot_LVLH = Rho_LVLH(4:6);

% Retrieve Target State in MCI
Xt_MCI = ppsval(XtPPsMCI, t);
COEt = rvPCI2COE(Xt_MCI, muM);
MEEt = COE2MEE(COEt);

% Convert Rho state into MCI
Rho_MCI = rhoLVLH2MCI(Rho_LVLH, Xt_MCI, t, EarthPPsMCI, SunPPsMCI, MoonPPsECI, muE, muS, deltaE, psiM, deltaM);

% Compute Chaser State in MCI
Xc_MCI = Xt_MCI + Rho_MCI;
COEc = rvPCI2COE(Xc_MCI, muM);
MEEc = COE2MEE(COEc);

% Target's Third, Fourth Body and Moon Harmonics Perturbing Accelerations
a34Bt = ThirdFourthBody(MEEt, t, EarthPPsMCI, SunPPsMCI, muE, muS);
aG_Mt = MoonHarmPerts(MEEt, MoonPPsECI, t, muM, deltaE, psiM, deltaM);
apt_LVLH = a34Bt + aG_Mt;

% Chaser's Third, Fourth Body and Moon Harmonics Perturbing Accelerations
a34Bc = ThirdFourthBody(MEEc, t, EarthPPsMCI, SunPPsMCI, muE, muS);
aG_Mc = MoonHarmPerts(MEEc, MoonPPsECI, t, muM, deltaE, psiM, deltaM);
apc_LVLH = a34Bc + aG_Mc;

% Convert Perturbating Accelerations into MCI
a_t = COEt(1);
e_t = COEt(2);
incl_t = COEt(3);
Omega_t = COEt(4);
omega_t = COEt(5);
nu_t = COEt(6);
theta_tt = COEt(5) + COEt(6);

R_MCI2LVLH = R3(theta_tt)*R1(incl_t)*R3(Omega_t);
apt_MCI = R_MCI2LVLH'*apt_LVLH;
apc_MCI = R_MCI2LVLH'*apc_LVLH;

% Compute Angular Velocity of LVLH wrt MCI = omega_LVLH
rt_MCI = Xt_MCI(1:3);
vt_MCI = Xt_MCI(4:6);

pt = a_t*(1-e^2);
rt = norm(rt_MCI);
ht_vect = cross(rt_MCI, vt_MCI);
ht = norm(ht_vect);

incl_t_dot = rt*apt_LVLH(3)/ht * cos(theta_tt);
Omega_t_dot = rt*apt_LVLH(3)/ht *sin(theta_tt)/sin(incl_t);
omega_t_dot = -rt*apt_LVLH(3)*sin(theta_tt)*cos(incl_t)/(ht*sin(incl_t)) - apt_LVLH(1)/e_t*cos(nu_t)*sqrt(pt/muM) + apt_LVLH(2)*sqrt(pt/muM)*sin(nu_t)*(e_t*cos(nu_t)+2)/(e_t*(1+e_t*cos(nu_t)));
nu_t_dot = sqrt(muM/pt^3)*(1+e_t*cos(nu_t))^2 + apt_LVLH(1)/e_t*cos(nu_t)*sqrt(pt/muM) - apt_LVLH(2)*sqrt(pt/muM)*sin(nu_t)*(e_t*cos(nu_t)+2)/e_t*(1+e_t*cos(nu_t));
theta_tt_dot = omega_t_dot + nu_t_dot;

omega_LVLH = [Omega_t_dot*sin(incl_t)*sin(theta_tt) + incl_t_dot*cos(theta_tt); ...
              Omega_t_dot*sin(incl_t)*cos(theta_tt) - incl_t_dot*sin(theta_tt); ...
              Omega_t_dot*cos(incl_t) + theta_tt_dot];


% Assign State Derivatives
dRho_LVLH(1:3) = rhodot_LVLH;
dRho_LVLH(4:6) = -2*cross(omega_LVLH, rhodot_LVLH) - cross(omegadot_LVLH, rho_LVLH) - cross(omega_LVLH, cross(omega_LVLH, rho_LVLH)) + ...
            muM/rt^3*((q*(2+q+(1+q)^(1/2)))/((1+q)^(3/2)*((1+q)^(1/2)+1)))*rt_MCI - muM/rc^3*rho_MCI + apc_MCI - apt_MCI;

end


