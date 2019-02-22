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

% Ln k0 dimensionless
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

zspan = 10:0.1:20;               % m
%solve odes
y0 = [0; 0; 0; 0; 0; 600; 130000];
[z, y] = ode45(@odefun, zspan, y0)

%% function space
function dydz = odefun(z, y)
global c
% arrhenius values for rates
k1 = exp((c.lnk0_1)-(c.ER_1/y(6)));
k2 = exp((c.lnk0_2)-(c.ER_2/y(6)));
k3 = exp((c.lnk0_3)-(c.ER_3/y(6)));
k4 = exp((c.lnk0_4)-(c.ER_4/y(6)));
k5 = exp((c.lnk0_5)-(c.ER_5/y(6)));

%initial moles in feed
c.n_oxi = 1.338;
c.n_o2i = 27.84;
c.n_n2i = 104.75;

% molar calculations
n_ox = c.n_oxi - y(1) - y(2) - y(3); % moles of oxylene 
n_o2 = c.n_o2i - 3*y(1) - 6.5*y(2) - 10.5*y(3) - 3.5*y(4) - 7.5*y(5); % moles of oxygen
n_pa = y(1) - y(4) - y(5); % moles of phthalic anhydride
n_w = 3*y(1) + 5*y(2) + 5*y(3)  - 2*y(4) + 2*y(5); % moles of water
n_co = 8*y(2) + 8*y(4); % moles of CO
n_co2 = 8*y(3) + 8*y(5); % moles of CO2

% total molar flowrate
nt = c.n_oxi + c.n_o2i + c.n_n2i + 5.5*y(2) + 1.5*y(3) + 5.5*y(4) + 1.5*y(5);

% volumetric flowrate
vt = (nt*8.314*y(6))/y(7);

% concentration calculations
C_ox = n_ox/vt;
C_o2 = n_o2/vt;
C_pa = n_pa/vt;
C_w = n_w/vt;
C_co = n_co/vt;
C_co2 = n_co2/vt;

% rates of reaction
r1 = (k1*power(C_o2, c.n)*c.b1*C_ox)/(1+c.b1*C_ox);
r2 = (k2*power(C_o2, c.n)*c.b2*C_ox)/(1+c.b2*C_ox);
r3 = (k3*power(C_o2, c.n)*c.b2*C_ox)/(1+c.b2*C_ox);
r4 = (k4*power(C_o2, c.n)*C_pa*c.b3);
r5 = (k5*power(C_o2, c.n)*C_pa*c.b3);

dydz = zeros(7,1); % initialising array of ODEs

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

% energy balance 
c.Tw = 610; % wall temp, K
c.a = 1.0310; 
c.b = -5*power(10,-5); 
c.c = 2.881*power(10,-7); 
c.d = -1.025*power(10,-10); % constants for temperature dependence of Cp
c.U = 0.096; %kW m^-2 K^-1
Q = c.A*c.U*(y(6)-c.Tw); % Q=A*U*(T-Tw), kW
cp = @(x) c.a + c.b*x + c.c*(power(x,2)) + c.d*(power(x,3));

% temperature (something in this
% equation may be causing the
% imaginary parts?
dydz(6) = (-Q-(r1*c.H1+r2*c.H2+r3*c.H3+r4*c.H4+r5*c.H5))*c.eps*c.A/(nt*cp(y(6)));

% pressure
dydz(7) = 1.3*power(10,5)-c.rho_c*(1-c.eps)*c.g*z;
end

