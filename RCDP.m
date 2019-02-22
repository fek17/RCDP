%% CE2-03-2 Group 6
clear, clc, close all

% setup for constants
global c

% define constants
c.b1 = 429;       % m^3.kmol^-1
c.b2 = 2024;      % m^3.kmol^-1
c.b3 = 1.0001;    % m^3.kmol^-1
c.n = 0.428;      % dimensionless
c.A = pi*power(0.025,2)/4;  % cross sectional area of reactor, m^2

% Ln k0
c.lnk0_1 = 19.84; 
c.lnk0_2 = 20.34;
c.lnk0_3 = 20.86;
c.lnk0_4 = 18.34;
c.lnk0_5 = 19.51;

% E/R (K)
c.ER_1 = 6920; 
c.ER_2 = 10436;
c.ER_3 = 10436;
c.ER_4 = 11457;
c.ER_5 = 11457;

% catalyst density
c.rho_c=1300;               % kg.m^-3
% bed voidage (epsilon)
c.eps=0.5;                  % dimensionless
% acceleration due to gravity
c.g=9.81;                   % m.s^-2
% initial pressure
c.Po=1.3*1.01*10^5;         % kg.m.s^-2; Po=Initial pressure

%% plot pressure vs riser height

% riser height
z = 0:0.1:10;               % m
plot(z,P(z))
xlabel('Riser Height (m)')
ylabel('Pressure (kg.m.s^{-2})')

%trying to solve odes
y0 = [0, 0, 0, 0, 0, 600];

%% function space

% pressure function
function p = P(z)
global c
p = c.Po - c.rho_c*(1-c.eps)*c.g*z;  % Pa
end

% differential
function dydz = odefun(z, y)
global c
% arrhenius values for rates
k1 = exp((c.lnk0_1)-(c.ER_1/T));
k2 = exp((c.lnk0_2)-(c.ER_2/T));
k3 = exp((c.lnk0_3)-(c.ER_3/T));
k4 = exp((c.lnk0_4)-(c.ER_4/T));
k5 = exp((c.lnk0_5)-(c.ER_5/T));

% rates of reaction
r1 = (k1*power(C_o2, c.n)*c.b1*C_ox)/(1+c.b1*C_ox);
r2 = (k2*power(C_o2, c.n)*c.b2*C_ox)/(1+c.b2*C_ox);
r3 = (k3*power(C_o2, c.n)*c.b2*C_ox)/(1+c.b2*C_ox);
r4 = (k4*power(C_o2, c.n)*C_pa*c.b3);
r5 = (k5*power(C_o2, c.n)*C_pa*c.b3);

dydz = zeros(6,1); % initialising array of ODEs

% dEi/dz = ri * epsilon * area
dydz(1) = r1*c.eps*c.A; 
dydz(2) = r2*c.eps*c.A;
dydz(3) = r3*c.eps*c.A;
dydz(4) = r4*c.eps*c.A;
dydz(5) = r5*c.eps*c.A;

% enthalpies of reaction (kJ/kmol)
c.H1 = -128300; 
c.H2 = -1526200;
c.H3 = -3273600;
c.H4 = -265600;
c.H5 = -1398000;

% energy balance (UNFINISHED!)
c.Tw = 610; % Twall temp, K
Q = c.A*0.096*(T-c.Tw); % Q=A*U*(T-Tw)
dydz(6) = Q-(r1*c.H1+r2*c.H2+r3*c.H3+r4*c.H4+r5*c.H5);
end

