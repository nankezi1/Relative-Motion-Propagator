function G = get_G(MEE, muM, eta)
% Description: this function defines the G matrix necessary in the
% propagation of Modified Equinoctial Elements.
% 
% Inputs:
% MEE = [x1, x2, x3, x4, x5, x6]
% 
% Outputs:
% G = gravitational model matrix

% Retrieve Data from MEE
x1 = MEE(1);
x2 = MEE(2);
x3 = MEE(3);
x4 = MEE(4);
x5 = MEE(5);
x6 = MEE(6);

G = sqrt(x1 / muM) * [0    2*x1/eta                               0;
    sin(x6)     ((1 + eta)*cos(x6) + x2) / eta    -(x4*sin(x6) - x5*cos(x6))*x3 / eta;
    -cos(x6)    ((1 + eta)*sin(x6) + x3) / eta    (x4*sin(x6) - x5*cos(x6))*x2 / eta;
       0                   0                      (1 + x4^2 + x5^2) / (2*eta)*cos(x6);
       0                   0                      (1 + x4^2 + x5^2) / (2*eta)*sin(x6);
       0                   0                      (x4*sin(x6) - x5*cos(x6)) / eta];

end