%% CE2-03-2 Group 6
clear, clc, close all

% setup for constants
global c

% define constants
c.b1 = 429;       % m^3.kmol^-1
c.b2 = 2024;      % m^3.kmol^-1
c.b3 = 1.0001;    % m^3.kmol^-1
c.n = 0.428;      % dimensionless

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

%% function space

% pressure function
function p = P(z)
global c
p = c.Po - c.rho_c*(1-c.eps)*c.g*z;  % Pa
end

% equations
function dydz = odefcn(z, y)
global c

k1 = exp((c.lnk0_1)-(c.ER_1/T));
k2 = exp((c.lnk0_2)-(c.ER_2/T));
k3 = exp((c.lnk0_3)-(c.ER_3/T));
k4 = exp((c.lnk0_4)-(c.ER_4/T));
k5 = exp((c.lnk0_5)-(c.ER_5/T));
dydz = zeros(6,1); % initialising array of extent of reaction
dydz(1) = (k1*power(C_o2, c.n)*c.b1*C_ox)/(1+c.b1*C_ox); % dE1/dz = r1
dydz(2) = (k2*power(C_o2, c.n)*c.b2*C_ox)/(1+c.b2*C_ox);
dydz(3) = (k3*power(C_o2, c.n)*c.b2*C_ox)/(1+c.b2*C_ox);
dydz(4) = (k4*power(C_o2, c.n)*C_pa*c.b3);
dydz(5) = (k5*power(C_o2, c.n)*C_pa*c.b3);
end

