%% CE2-03-2 Group 6
%  reactor

format long

% constants, variables, temporary & output things
global c v t o

% bring in constants.m only if sensitivity test flag DNE or is false
if ~isfield(c,'Canary') || c.Canary == false
    constants
end

% cross sectional area of reactor
c.A = pi*c.Dia^2/4;   % m^2

%% feed

% mole fractions
c.S.x0 = zeros(7,1);
c.S.('x0')('OX') = c.fOX;
c.S.('x0')('N2') = (1-c.fOX)*0.79;
c.S.('x0')('O2') = (1-c.fOX)*0.21;

% bring in variable
c.f.massFlux = c.fMFlux;

% total Mw of feed
c.f.Mw = dot(c.S.x0,c.S.Mw);            % kg.kmol^{-1}

% feed mass flow
c.f.massFlow = c.f.massFlux * c.A;      % kg.h^-1
% *assumption klaxon* ? superficial mass flux so use total area of column

% feed molar flow
c.f.n0 = c.f.massFlow / c.f.Mw ;        % kmol.h^{-1}

% component molar flows
c.S.n0 = c.f.n0 * c.S.x0;
c.S.Properties.VariableUnits{'n0'} = 'kmol.h^{-1}';

%% solver

% height of reactor
t.zspan = 0:0.01:5;           % m

% initial conditions
t.y0 = [0; 0; 0; 0; 0; c.T0; c.P0]; 


% solve odes
[t.z, t.y] = ode15s(@odefun, t.zspan, t.y0);

% dump output in table
v = table;
v(:,1:7) = array2table(t.y);
v.Properties.VariableNames = {'xi_1'        'xi_2'        'xi_3'        'xi_4'        'xi_5'        'T' 'P' };
v.Properties.VariableUnits = {'kmol.h^{-1}' 'kmol.h^{-1}' 'kmol.h^{-1}' 'kmol.h^{-1}' 'kmol.h^{-1}' 'K' 'Pa'};
v.z = t.z;

%% molar flows

% for each z position in output
for i = 1:size(v,1)
    
    % do material balances
    t.S = reactorMB(c.S, v{i,{'xi_1' 'xi_2' 'xi_3' 'xi_4' 'xi_5'}}, v.T(i), v.P(i));

    % store each calculated flow
    for j = 1:numel(c.species)
        n_j = sprintf('n_%s',c.species{j});
        
        % initialise columns in v on first run
        if i == 1
            v.(n_j) = zeros(numel(t.zspan),1);
        end
        
        % store output in table
        v.(n_j)(i) = t.S{c.species{j},'n'};
    end
    
end
clear i j n_j

%% KPIs

% yield of PA wrt. OX
v.Y_PA_OX = v.n_PA / v.n_OX(1);

% OX conversion
v.f_OX = (v.n_OX(1) - v.n_OX)/v.n_OX(1);

% selectivity of PA (wrt. OX *assumption*)
v.S_PA_OX = v.n_PA ./ (v.n_OX(1) - v.n_OX);

% selectivity of CO (wrt. OX *assumption*)
v.S_CO_OX = v.n_CO ./ (v.n_OX(1) - v.n_OX);

% CO/CO2 fraction
v.CO_CO2 = v.n_CO ./ v.n_CO2;

% output (for sensitivity analysis)
for i = 1:numel(KPI)

    % max value & z position
    [o.(KPI{i}).max.val, t.I_z] = max(v.(KPI{i}));
    o.(KPI{i}).max.z = v.z(t.I_z);
    
end
clear i

%% function space

% ode function for solver
function dydz = odefun(z, y)

% bring in constants
global c
xi = y(1:5);
T = y(6);
P = y(7);

% kinetic constant
k = exp( c.RX.lnk0 - c.RX.ER/T ); % h^{-1}.kmol^{1-n}.m^{3n}.kg_cat^{-1}

% do material balances
S = reactorMB(c.S, xi, T, P);

% rates of reaction ( kmol.h^{-1}.kg_cat^{-1} )
r_c(1)   = ( k(1)   * S.('C')('O2')^c.n * c.b1 * S.('C')('OX') )/( 1 + c.b1 * S.('C')('OX') );
r_c(2:3) = ( k(2:3) * S.('C')('O2')^c.n * c.b2 * S.('C')('OX') )/( 1 + c.b2 * S.('C')('OX') );
r_c(4:5) = ( k(4:5) * S.('C')('O2')^c.n * S.('C')('PA') * c.b3 );

% rates of reaction ( kmol.h^{-1}.m^{-3} )
r = r_c * c.rho_c * (1-c.eps)/c.eps;

% extents [ kmol.h^{-1}.m^{-1} ]
dxidz = r * c.eps * c.A;

% heat capacity function
cp = @(T) c.a + c.b*T + c.c*T^2 + c.d*T^3; % kJ kg^-1 K^-1

% temperature
dTdz = ( - c.U * pi * c.Dia * (T - c.Tw) - dot(c.RX.h,dxidz) )/( cp(T) * c.f.massFlow );

% pressure
dPdz = - c.rho_c * (1-c.eps) * c.g;

% output derivatives
dydz = [ dxidz dTdz dPdz ]';

end

% material balances
function [S] = reactorMB(S, xi, T, P)

% inputs
%   S   (table)
%       rows: components j
%       cols: properties
%           n0  initial molar flow
%                   kmol.h^{-1}
%   xi  (vector) extents of reaction i (1-5)
%           kmol.h^{-1}
%   T   temperature
%           K
%   P   pressure
%           Pa
%
% outputs
%   S   (table)
%       rows: components j
%       cols: properties
%           n0  initial molar flow
%                   kmol.h^{-1}
%           n   molar flow
%                   kmol.h^{-1}
%           C   concentration
%                   kmol.m^{-3}

global c

% individual molar flowrate
S.n = zeros(7,1);
S.Properties.VariableUnits{'n'} = 'kmol.h^{-1}';
S.('n')('OX')  = S.('n0')('OX') - xi(1) - xi(2) - xi(3);
S.('n')('O2')  = S.('n0')('O2') - 3*xi(1) - 6.5*xi(2) - 10.5*xi(3) - 3.5*xi(4) - 7.5*xi(5);
S.('n')('PA')  = S.('n0')('PA') + xi(1) - xi(4) - xi(5);
S.('n')('H2O') = S.('n0')('H2O') + 3*xi(1) + 5*xi(2) + 5*xi(3) + 2*xi(4) + 2*xi(5);
S.('n')('CO')  = S.('n0')('CO') + 8*xi(2) + 8*xi(4);
S.('n')('CO2') = S.('n0')('CO2') + 8*xi(3) + 8*xi(5);
S.('n')('N2')  = S.('n0')('N2');

% total volumetric flowrate
vt = (sum(S.n) * c.R * T)/P * 10^3; % m^3.h^{-1}

% concentration calculations
S.C = S.n ./ vt;
S.Properties.VariableUnits{'C'} = 'kmol.m^{-3}';

end