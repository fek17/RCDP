%% CE2-03-2 Group 6
%  constants

clear global c
global c

%% general

% diameter of reactor 
c.D = 0.025;        % m
% cross sectional area of reactor
c.A = pi*c.D^2/4;   % m^2
% acceleration due to gravity
c.g = 9.81;         % m.s^-2
% initial pressure
c.P0 = 1.3*10^5;    % Pa
% initial temperature
c.T0 = 600;         % K
% ideal gas constant
c.R = 8.3144598;    % J.K^{-1}.mol^{-1}

%% kinetics

% constants in rate expressions
c.b1 = 429;         % m^3.kmol^-1
c.b2 = 2024;        % m^3.kmol^-1
c.b3 = 1.0001;      % m^3.kmol^-1
c.n = 0.428;        % dimensionless

% catalyst density
c.rho_c = 1300;     % kg.m^-3
% bed voidage (epsilon)
c.eps = 0.5;        % dimensionless

% table of reaction data ( ln k_0  E/R ?h )
c.RX = table;
c.RX.lnk0 = [19.24 20.04 20.36 28.55 28.45]';
c.RX.ER = [14987 15862 15862 16040 16040]';
c.RX.h = [-1283000 -1526200 -3273600 -265600 -1398000]';
c.RX.Properties.VariableUnits = {'dimensionless','K','kJ.kmol^{-1}'};

%% species

% setup table
c.species = {'OX' 'PA' 'CO2' 'CO' 'H2O' 'O2' 'N2'};
c.S = table('RowNames',c.species);

% molecular weights
c.S.Mw = [106.1602 148.1100 44.0095 28.0101 18.0153 31.9988 28.0134]';
c.S.Properties.VariableUnits{'Mw'} = 'kg.kmol^{-1}';

%% feed

% feed mole fractions
c.S.x0 = zeros(7,1);
c.S{'OX','x0'} = 0.01;
c.S{'N2','x0'} = 0.78;
c.S{'O2','x0'} = 0.21;

% feed mass flux
c.f.massFlux = 2500;                    % kg.m^-2.h^-1

%% energy

% heat capacity constants
c.a = 1.0310; 
c.b = -5e-5; 
c.c = 2.881e-7; 
c.d = -1.025e-10;

% wall temp
c.Tw = 610; % K

% wall heat transfer coefficient
c.U = 0.096*60^2; %kJ h-1 m^-2 K^-1