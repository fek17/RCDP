%% CE2-03-2 Group 6

% define constants
b1 = 429;       % m^3.kmol^-1
b2 = 2024;      % m^3.kmol^-1
b3 = 1.0001;    % m^3.kmol^-1
n = 0.428;      % dimensionless

% equations
y=2;

%plot pressure vs riser height

rho_c=1300;              % kg.m^-3; rho_c=catalyst density
eps=0.5;                 % dimensionless; eps=bed voidage
g=9.81;                  % m.s^-2; g=gravitational acceleration constant
z=0:0.1:10;              % z=riser height
Po=1.3*1.01*10^5         %kg.m.s^-2; Po=Initial pressure
P=Po-rho_c*(1-eps)*g*z;  %P=pressure
plot(z,P)
xlabel('Riser Height (m)')
ylabel('Pressure (kg.m.s^{-2})')


