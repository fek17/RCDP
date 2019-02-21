%% CE2-03-2 Group 6

% define constants
b1 = 429;       % m^3.kmol^-1
b2 = 2024;      % m^3.kmol^-1
b3 = 1.0001;    % m^3.kmol^-1
n = 0.428;      % dimensionless
lnk0_1 = 19.84; 
lnk0_2 = 20.34;
lnk0_3 = 20.86;
lnk0_4 = 18.34;
lnk0_5 = 19.51;
ER_1 = 

% plot pressure vs riser height



% equations
function dEdz = odefcn(P, T)

k1 =
k2 =
k3 =
k4 =
k5 =
dEdz = zeros(5,1); % initialising array of extent of reaction
dEdz(1) = (k1*power(C_o2, n)*b1*C_ox)/(1+b1*C_ox); % dE1/dz = r1
dEdz(2) = (k2*power(C_o2, n)*b2*C_ox)/(1+b2*C_ox);
dEdz(3) = (k3*power(C_o2, n)*b2*C_ox)/(1_b2*C_ox);
dEdz(4) = (k4*power(C_o2, n)*C_pa*b3)
dEdz(5) = (k5*power(C_o2, n)*C_pa*b3)
end