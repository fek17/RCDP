%% CE2-03-2 Group 6
%  Plots

global v

% run reactor if we don't have the data already
if ~istable(v)
    reactor
end

% export?
global printFlag
printFlag = false;

%% extents of reaction
figure
hold on
for i = 1:size(c.RX,1)
    t.this = sprintf('xi_%u',i);
    plot(v.z,v.(t.this),'DisplayName',['\' t.this])
end
ylabel('\xi_i / kmol.h^{-1}');
xlabel('z / m');
legend('Location','southoutside','Orientation','horizontal','NumColumns',3);
legend('boxoff');

figExport(8,10,'overview-extents');

%% molar flows (except inert N2)

figure
hold on
for i = 1:numel(c.species)
    % remove inert
    if strcmp(c.species{i},'N2') == false
        plot(v.z,v.(['n_' c.species{i}]),'DisplayName',['n_{' c.species{i} '}'])
    end
end
ylabel('n_j / kmol.h^{-1}');
xlabel('z / m');
ylim([1e-6 1e-2]);
legend('Location','southoutside','Orientation','horizontal','NumColumns',3);
legend('boxoff');

figExport(8,10,'overview-molar-flows');

% log scale version
set(gca, 'YScale', 'log');
figExport(8,10,'overview-molar-flows-log');

%% yield of PA wrt OX

figure
plot(v.z,v.Y_PA_OX)
ylabel('Y_{PA,OX}^{OV}');
xlabel('z / m');

figExport(8,8,'yield-PA-OX');

%% temperature

figure
plot(v.z,v.T)
ylabel('T / K');
xlabel('z / m');

figExport(8,8,'temp');

%% pressure

figure
plot(v.z,v.P)
ylabel('P / Pa');
xlabel('z / m');

figExport(8,8,'pressure');

%% combined T P

figure
yyaxis left
plot(v.z,v.T)
ylabel('T / K');
xlabel('z / m');
yyaxis right
plot(v.z,v.P)
ylabel('P / Pa');

figExport(8,8,'combined-temp-pressure');

%% OX conversion

figure
plot(v.z,v.f_OX)
ylabel('f_{OX}^{OV}');
xlabel('z / m');

figExport(8,8,'conversion-OX');

%% combined Y f

figure
yyaxis left
plot(v.z,v.Y_PA_OX)
ylabel('Y_{PA,OX}^{OV}');
xlabel('z / m');
yyaxis right
plot(v.z,v.f_OX)
ylabel('f_{OX}^{OV}');

figExport(8,8,'combined-yield-PA-OX-conversion-OX');

%% selectivity of (PA, CO) wrt. OX

figure
yyaxis left
plot(v.z,v.S_PA_OX)
ylabel('S_{PA,OX}');
xlabel('z / m');
yyaxis right
plot(v.z,v.S_CO_OX)
ylabel('S_{CO,OX}');

figExport(8,8,'selectivity-PA-CO-OX');

%% ratio CO/CO2

figure
plot(v.z,v.CO_CO2)
ylabel('n_{CO}/n_{CO2}');
xlabel('z / m');

figExport(8,8,'ratio-CO-CO2');
