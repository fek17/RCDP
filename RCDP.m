%% CE2-03-2 Group 6
clear, clc, close all
format long

%% constants
global c

% diameter of reactor 
c.D = 0.025; % m
% cross sectional area of reactor
c.A = pi*power(c.D,2)/4;  % m^2
% acceleration due to gravity
c.g = 9.81;         % m.s^-2
% initial pressure
c.Po = 1.3*10^5;    % Pa
% ideal gas constant
c.R = 8.3144598;    % J.K^{-1}.mol^{-1}

%% kinetics

% define constants
c.b1 = 429;         % m^3.kmol^-1
c.b2 = 2024;        % m^3.kmol^-1
c.b3 = 1.0001;      % m^3.kmol^-1
c.n = 0.428;        % dimensionless

% catalyst density
c.rho_c = 1300;     % kg.m^-3
% bed voidage (epsilon)
c.eps = 0.5;        % dimensionless

% Ln k0 dimensionless
c.lnk0_1 = 19.24; 
c.lnk0_2 = 20.04;
c.lnk0_3 = 20.36;
c.lnk0_4 = 28.55;
c.lnk0_5 = 28.45;

% E/R (K)
c.ER_1 = 14987; 
c.ER_2 = 15862;
c.ER_3 = 15862;
c.ER_4 = 16040;
c.ER_5 = 16040;

%% feed

% feed species
c.f.T = table('RowNames',{'OX';'N2';'O2'});
% mole fractions
c.f.T.x_i = [0.01; 0.78; 0.21];
% molecular weight
c.f.T.Mw = [106.1602; 28.0134; 31.9988];
c.f.T.Properties.VariableUnits{'Mw'} = 'kg.kmol^{-1}';

% total Mw of feed
c.f.Mw = dot(c.f.T.x_i,c.f.T.Mw);       % kg.kmol^{-1}

% feed mass flux
c.f.massFlux = 2500;                    % kg.m^-2.h^-1

% feed mass flow
c.f.massFlow = c.f.massFlux * c.A;      % kg.h^-1
% *assumption klaxon* ? superficial mass flux so use total area of column

% feed molar flow
c.f.Mf = c.f.massFlow / c.f.Mw ;   % kmol.h^{-1}

% component molar flows
c.f.T.Mf = c.f.Mf * c.f.T.x_i;
c.f.T.Properties.VariableUnits{'Mf'} = 'kmol.h^{-1}';

% extract initial feed values for solver ( kmol.h^{-1} )
c.n_oxi = c.f.T{'OX','Mf'};
c.n_o2i = c.f.T{'O2','Mf'};
c.n_n2i = c.f.T{'N2','Mf'}; % nitrogen in air, inert
c.n_pai = 0;
c.n_wi = 0;
c.n_coi = 0;
c.n_co2i = 0;

%% solver

% height of reactor
zspan = 0:0.1:20;           % m

% initial conditions
y0 = [0; 0; 0; 0; 0; 600; c.Po]; % 5x extents [kmol.h^{-1}], temperature [K], pressure [Pa]

% solve odes
[z, y] = ode45(@odefun, zspan, y0);

%% plot

figure
yyaxis left
hold on
for i = 1:5
    plot(z,y(:,i),'DisplayName',sprintf('\\xi_%u',i))
end
ylabel('\xi_i');
yyaxis right
plot(z,y(:,6),'DisplayName','T / K')
ylabel('T / K');
legend('Location','east');

%% function space

function dydz = odefun(z, y)

% bring in constants
global c

% kinetic constant ( h^{-1}.kmol^{1-n}.m^{3n}.kg_cat^{-1} )
k1 = exp((c.lnk0_1)-(c.ER_1/y(6)));
k2 = exp((c.lnk0_2)-(c.ER_2/y(6)));
k3 = exp((c.lnk0_3)-(c.ER_3/y(6)));
k4 = exp((c.lnk0_4)-(c.ER_4/y(6)));
k5 = exp((c.lnk0_5)-(c.ER_5/y(6)));

% molar flow calculations [ kmol.h^{-1} ]
n_ox = c.n_oxi - y(1) - y(2) - y(3); % moles of o-xylene 
n_o2 = c.n_o2i - 3*y(1) - 6.5*y(2) - 10.5*y(3) - 3.5*y(4) - 7.5*y(5); % moles of oxygen
n_pa = c.n_pai + y(1) - y(4) - y(5); % moles of phthalic anhydride
n_w = c.n_wi + 3*y(1) + 5*y(2) + 5*y(3) + 2*y(4) + 2*y(5); % moles of water
n_co = c.n_coi + 8*y(2) + 8*y(4); % moles of CO
n_co2 = c.n_co2i + 8*y(3) + 8*y(5); % moles of CO2
n_n2 = c.n_n2i;

% total molar flowrate [ kmol.h^{-1} ]
nt = c.n_oxi + c.n_o2i + c.n_n2i + 5.5*y(2) + 1.5*y(3) + 5.5*y(4) + 1.5*y(5);
% nt = n_ox + n_o2 + n_pa + n_w + n_co + n_co2 + n_n2; % alt. approach

% volumetric flowrate
vt = (nt * (c.R*y(6))/y(7)) * 10^3; % m^3.h^{-1}

% concentration calculations
C_ox = n_ox/vt;
C_o2 = n_o2/vt;
C_pa = n_pa/vt;
C_w = n_w/vt;
C_co = n_co/vt;
C_co2 = n_co2/vt;

% rates of reaction ( kmol.h^{-1}.kg_cat^{-1} )
r1 = (k1*power(C_o2, c.n)*c.b1*C_ox)/(1+c.b1*C_ox);
r2 = (k2*power(C_o2, c.n)*c.b2*C_ox)/(1+c.b2*C_ox);
r3 = (k3*power(C_o2, c.n)*c.b2*C_ox)/(1+c.b2*C_ox);
r4 = (k4*power(C_o2, c.n)*C_pa*c.b3);
r5 = (k5*power(C_o2, c.n)*C_pa*c.b3);

% initialising array for derivatives
dydz = zeros(7,1);

% d(kg_c)/dz
c.dkgcdz = c.A * (1-c.eps) * c.rho_c;

% void area
c.A_v = c.eps * c.A;

% extents [ kmol.h^{-1}.m^{-1} ]
dydz(1) = r1 * c.A * c.rho_c * (1-c.eps);
dydz(2) = r2 * c.A * c.rho_c * (1-c.eps);
dydz(3) = r3 * c.A * c.rho_c * (1-c.eps);
dydz(4) = r4 * c.A * c.rho_c * (1-c.eps);
dydz(5) = r5 * c.A * c.rho_c * (1-c.eps);

% enthalpies of reaction (kJ/kmol)
c.H1 = -1283000; 
c.H2 = -1526200;
c.H3 = -3273600;
c.H4 = -265600;
c.H5 = -1398000;

% heat capacity constants
c.a = 1.0310; 
c.b = -5e-5; 
c.c = 2.881e-7; 
c.d = -1.025e-10;

% heat capacity function
cp = @(T) c.a + c.b*T + c.c*(power(T,2)) + c.d*(power(T,3)); % kJ kg^-1 K^-1

% energy balance 

c.Tw = 610; % wall temp, K

c.U = 0.096*3600; %kJ h-1 m^-2 K^-1
Q = c.D*pi*c.U*(y(6)-c.Tw); % Q=A*U*(T-Tw), kJ h-1

% temperature
dydz(6) = ((-Q-(y(1)*c.H1+y(2)*c.H2+y(3)*c.H3+y(4)*c.H4+y(5)*c.H5))*c.eps*c.A)/(c.mt*cp(y(6)));

% pressure
dydz(7) = 1.3*power(10,5)-(c.rho_c*(1-c.eps)*c.g*z);
end

