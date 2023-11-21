function dRHO_LVLH = DynamicalModelRM(t, RHO_LVLH, EarthPPsMCI, SunPPsMCI, muE, ...
                                      muS, tspan, MoonPPsECI, deltaE, psiM, deltaM, ...
                                      t0, tf, XtPPsMCI, omegaPPsLVLH)

% Description: this is the function with the Dynamical Model for the
% Relative Motion.
% 
% Inputs:
% t = epoch
% RHO_LVLH = [rho, rhodot] state in LVLH
% EarthPPsMCI = pp struct with the interpolation of the Earth State in MCI
% SunPPsMCI = pp struct with the interpolation of the Sun State in MCI
% muE = Earth's gravitational parameter in canonical units
% muS = Sun's gravitational parameter in canonical units
% tspan = time span
% MoonPPsECI = pp struct with the interpolation of the Moon State in ECI
% deltaE = Earth's Ecliptic Obliquity
% psiM = Moon's Precession Angle
% deltaM = Moon's Equator Obliquity
% t0 = initial time
% tf = final time
% XtPPsMCI = pp struct for the interpolation of target state in MCI
% omegaPPsLVLH = pp struct for the interpolation of LVLH frame angular velocity
% 
% Outputs:
% dRHO_LVLH = derivatives [rho, rhodot] state in LVLH

% Retrieve Global Variables
global muM Rm DU TU pbar log

% Initialize Derivatives Vector
dRHO_LVLH = zeros(6, 1);

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
rho_LVLH = RHO_LVLH(1:3);
rhodot_LVLH = RHO_LVLH(4:6);

% Retrieve Target State in MCI
Xt_MCI = ppsval(XtPPsMCI, t);
COEt = rvPCI2COE(Xt_MCI', muM)';
MEEt = COE2MEE(COEt')';
rt_MCI = Xt_MCI(1:3);
rt = norm(rt_MCI);

% Convert RHO state into MCI
RHO_MCI = rhoLVLH2MCI(RHO_LVLH, Xt_MCI, t, EarthPPsMCI, SunPPsMCI, MoonPPsECI, muE, muS, deltaE, psiM, deltaM);
rho_MCI = RHO_MCI(1:3);

% Compute Chaser State in MCI
Xc_MCI = Xt_MCI + RHO_MCI;
COEc = rvPCI2COE(Xc_MCI', muM)';
MEEc = COE2MEE(COEc')';
rc_MCI = Xc_MCI(1:3);
rc = norm(rc_MCI);

q = (dot(rho_MCI, rho_MCI) + 2*dot(rho_MCI, rt_MCI))/rt^2;

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
[omega_r, omegadot_r] = PolyEval(t, tspan, flip(omegaPPsLVLH(1).coefs, 2));
[omega_t, omegadot_t] = PolyEval(t, tspan, flip(omegaPPsLVLH(2).coefs, 2));
[omega_h, omegadot_h] = PolyEval(t, tspan, flip(omegaPPsLVLH(3).coefs, 2));
omega_LVLH = [omega_r, omega_t, omega_h]';
omegadot_LVLH = [omegadot_r, omegadot_t, omegadot_h]';

% Assign State Derivatives
dRHO_LVLH(1:3) = rhodot_LVLH;
dRHO_LVLH(4:6) = -2*cross(omega_LVLH, rhodot_LVLH) - cross(omegadot_LVLH, rho_LVLH) - cross(omega_LVLH, cross(omega_LVLH, rho_LVLH)) + ...
            muM/rt^3*((q*(2+q+(1+q)^(1/2)))/((1+q)^(3/2)*((1+q)^(1/2)+1)))*rt_MCI - muM/rc^3*rho_MCI + apc_MCI - apt_MCI;

upd = dRHO_LVLH./RHO_LVLH;
tol = 1;

for i = 1 : length(upd)
    if abs(upd(i)) > tol
        fprintf(log, 'Update to State Ratio:\n[%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f]\n', upd);
        fprintf(log, 'Time Elapsed: %02d days, %02d hrs, %02d mins\n\n', tDAY, tHR, tMIN);
        break
    end
end

end


