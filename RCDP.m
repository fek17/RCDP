%% CE2-03-2 Group 6

% define constants
b1 = 429;       % m^3.kmol^-1
b2 = 2024;      % m^3.kmol^-1
b3 = 1.0001;    % m^3.kmol^-1
n = 0.428;      % dimensionless

%Ln k0
lnk0_1 = 19.84; 
lnk0_2 = 20.34;
lnk0_3 = 20.86;
lnk0_4 = 18.34;
lnk0_5 = 19.51;

%E/R (K)
ER_1 = 6920; 
ER_2 = 10436;
ER_3 = 10436;
ER_4 = 11457;
ER_5 = 11457;

%plot pressure vs riser height
rho_c=1300;              % kg.m^-3; rho_c=catalyst density
eps=0.5;                 % dimensionless; eps=bed voidage
g=9.81;                  % m.s^-2; g=gravitational acceleration constant
z=0:0.1:10;              % z=riser height (m)
Po=1.3*1.01*10^5         %kg.m.s^-2; Po=Initial pressure
P=Po-rho_c*(1-eps)*g*z;  %P=pressure
plot(z,P)
xlabel('Riser Height (m)')
ylabel('Pressure (kg.m.s^{-2})')

% equations
function dEdz = odefcn(P, T)

k1 = exp((lnk0_1)-(ER_1/T))
k2 = exp((lnk0_2)-(ER_2/T))
k3 = exp((lnk0_3)-(ER_3/T))
k4 = exp((lnk0_4)-(ER_4/T))
k5 = exp((lnk0_5)-(ER_5/T))
dEdz = zeros(5,1); % initialising array of extent of reaction
dEdz(1) = (k1*power(C_o2, n)*b1*C_ox)/(1+b1*C_ox); % dE1/dz = r1
dEdz(2) = (k2*power(C_o2, n)*b2*C_ox)/(1+b2*C_ox);
dEdz(3) = (k3*power(C_o2, n)*b2*C_ox)/(1_b2*C_ox);
dEdz(4) = (k4*power(C_o2, n)*C_pa*b3)
dEdz(5) = (k5*power(C_o2, n)*C_pa*b3)
end

