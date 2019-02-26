%% RCDP Matlab refresher
%  Andreas Richardson
clear, clc, close all

%% constants

V = 1;          % m^3
E_a = 41840;    % kJ.kmol^-1
k_0 = 148.4;    % (kmol.m^-3)^-2 .K^-1
DH_r = 104670;  % kJ.kmol^-1
C_p = 2346.5;   % kJ.M^-3.K^-1
R = 8.314;      % kJ.kmol^-1.K^-1

%% calculations

% rate constant (T in K)
k = @(T) k_0 * exp( - E_a / (R * T) ); % (kmol.m^-3)^-2 .K^-1

% initial values (C_A / kmol.m^-3 , T / K)
y_0 = [10.2; 673.15];

% time range
t_range = [0 1]; % s

% define anonymous function (only 2 inputs) to pass through constants
thisReactor = @(t, y) reactorDerivs(t, y, k, DH_r, C_p);

% solve away!
[t, y] = ode45(thisReactor, t_range, y_0);

%% plot

figure
yyaxis left
plot(t,y(:,1));
ylabel("C_A / kmol.m^{-3}");
yyaxis right
plot(t,y(:,2));
ylabel("T / K");

%% function space

function [dy] = reactorDerivs(t, y, k, DH_r, C_p)

% extract values from input
C_A = y(1); % kmol.m^-3
T = y(2);   % K

% calculate derivs
dC_A = -k(T) .* C_A.^3;   % kmol.m^-3.s^-1
dT = -DH_r./C_p .* dC_A;  % K.s^-1

% output
dy = [dC_A; dT];

end