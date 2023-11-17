function Rho_MCI = rhoLVLH2MCI(Rho_LVLH, Xt_MCI, t, EarthPPsMCI, SunPPsMCI, MoonPPsECI, muE, muS, deltaE, psiM, deltaM)
% Description: this function converts vectors in the MCI reference frame
% into the LVLH rotating reference frame centered on the Target trajectory.

% Recall Global Variables
global muM Rm DU TU

% Retrieve Data from Input
rho_LVLH = Rho_LVLH(1:3);
rhodot_LVLH = Rho_LVLH(4:6);

% Compute MCI2LVLH Rotation Matrix and its Derivative
[R, Rdot] = get_R_Rdot(Xt_MCI, t, EarthPPsMCI, SunPPsMCI, MoonPPsECI, muE, muS, deltaE, psiM, deltaM);

% Rotate both Position and Velocity
rho_MCI = R'*rho_LVLH;
rhodot_MCI = R'*(rhodot_LVLH - Rdot*R'*rho_LVLH);

% Assemble the MCI State
Rho_MCI = [rho_MCI; rhodot_MCI];

end