%% CE2-03-2 Group 6
%  plots for real-time controller test

close all

%% inlet temperature profiles

figure
plotStackedDataSeries(logsout,'T_r_in','T_c_in')

figExport(8,6,'temp-in')

%% PA flowrate

figure
plotDataSeries(logsout,'Outlet - PA flowrate');
plotBounds([1.5 1])

figExport(8,6,'PA-flow')

%% cooling medium outlet

figure
plotDataSeries(logsout,'T_c_out');
plotBounds([730 710 750])

figExport(8,6,'temp-out')

%% TC performance

figure
plotAdjacentDataSeries(logsout,'e_TC','op_TC')

figExport(16,4,'TC-performance')

%% FC performance

figure
plotAdjacentDataSeries(logsout,'e_FC','op_FC')

figExport(16,4,'FC-performance')

%% Smith predictor

figure
plotStackedDataSeries(logsout,'SP_out','T_c_out')

figExport(8,6,'sp')

%% plotter
function plotDataSeries(logsout,name)
series = getElementNames(logsout);
index = dot(1:numel(series),strcmp(name,series));
plot(logsout{index}.Values);
title('');
xlim([0 1010]);
end

function plotStackedDataSeries(logsout,name1,name2)
subplot(2,1,1)
plotDataSeries(logsout,name1)
plainTextY;
xlabel('')
subplot(2,1,2)
plotDataSeries(logsout,name2)
plainTextY;
end

function plotAdjacentDataSeries(logsout,name1,name2)
yyaxis left
plotDataSeries(logsout,name1)
plainTextY;
xlabel('')
yyaxis right
plotDataSeries(logsout,name2)
plainTextY;
end

function plotBounds(val)
x = [0 1010];
for i=1:numel(val)
    thisLine = line( [x(1); x(2)], val(i)*ones(2,1));
    if i == 1
        thisLine.Color = 'g';
        thisLine.LineStyle = '-.';
    else
        thisLine.Color = 'r';
        thisLine.LineStyle = '--';
    end
end
end

function plainTextY
ax = gca;
set(ax.YLabel,'Interpreter','none');
end

