function dRho_LVLH = DynamicalModelRM(t, Rho_LVLH, EarthPPsMCI, SunPPsMCI, muE, ...
                     muS, MoonPPsECI, deltaE, psiM, deltaM, t0, tf, XtPPsMCI, omegaPPsLVLH)
% Description: this is the function with the Dynamical Model for the
% Relavtive Motion.

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
[R_MCI2LVLH, ~] = get_R_Rdot(Xt_MCI, t, EarthPPsMCI, SunPPsMCI, MoonPPsECI, muE, muS, deltaE, psiM, deltaM);
apt_MCI = R_MCI2LVLH'*apt_LVLH;
apc_MCI = R_MCI2LVLH'*apc_LVLH;

% Compute Angular Velocity of LVLH wrt MCI and its derivative
omega_LVLH = ppsval(omegaPPsLVLH, t);



% Assign State Derivatives
dRho_LVLH(1:3) = rhodot_LVLH;
dRho_LVLH(4:6) = -2*cross(omega_LVLH, rhodot_LVLH) - cross(omegadot_LVLH, rho_LVLH) - cross(omega_LVLH, cross(omega_LVLH, rho_LVLH)) + ...
            muM/rt^3*((q*(2+q+(1+q)^(1/2)))/((1+q)^(3/2)*((1+q)^(1/2)+1)))*rt_MCI - muM/rc^3*rho_MCI + apc_MCI - apt_MCI;

end


