function RHO_LVLH = rhoMCI2LVLH(RHO_MCI, Xt_MCI, t, EarthPPsMCI, SunPPsMCI, MoonPPsECI, muE, muS, deltaE, psiM, deltaM)
% Description: this function converts vectors in the MCI reference frame
% into the LVLH rotating reference frame centered on the Target trajectory.

% Recall Global Variables
global muM Rm DU TU

% Retrieve Data from Input
rho_MCI = RHO_MCI(1:3);
rhodot_MCI = RHO_MCI(4:6);

% Compute MCI2LVLH Rotation Matrix and its Derivative
[R, Rdot] = get_R_Rdot(Xt_MCI, t, EarthPPsMCI, SunPPsMCI, MoonPPsECI, muE, muS, deltaE, psiM, deltaM);

% Rotate both Position and Velocity
rho_LVLH = R*rho_MCI;
rhodot_LVLH = R*rhodot_MCI + Rdot*rho_MCI;

% Assemble the LVLH State
RHO_LVLH = [rho_LVLH; rhodot_LVLH];

end