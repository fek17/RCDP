%% CE2-03-2 Group 6
%  sensitivity plots

% load data exported from sensitivity.m
if ~isstruct(s)
    load('s.mat')
end

% export?
global printFlag
printFlag = true;

%% b1 b2 b3 on Y_PA

i_kpi = 1;
fields = {'b1' 'b2' 'b3'};

multiPlot(s,fields,i_kpi,KPI)

figExport(15,4,'sensitivity/overview-b1-b2-b3')

%% ER 4/5 on Y_PA

i_kpi = 1;
fields = {'ER_4' 'ER_5'}; 

multiPlot(s,fields,i_kpi,KPI)

figExport(10,4,'sensitivity/overview-ER_4-ER_5')

%% dia on f_OX

i_kpi = 2;
fields = {'Dia'};

multiPlot(s,fields,i_kpi,KPI)

figExport(5,4,'sensitivity/overview-Dia')

%% eps on f_OX

i_kpi = 2;
fields = {'eps'}; 

multiPlot(s,fields,i_kpi,KPI)

figExport(5,4,'sensitivity/overview-eps')

%% fOX

i_kpi = 3;
fields = {'fOX'};

multiPlot(s,fields,i_kpi,KPI)

figExport(5,4,'sensitivity/overview-fOX')

%% P0

i_kpi = 1;
fields = {'P0'};

multiPlot(s,fields,i_kpi,KPI)

figExport(5,4,'sensitivity/overview-P0')

%% U

i_kpi = 2;
fields = {'U'};

multiPlot(s,fields,i_kpi,KPI)

yyaxis left
ylim([0.9 1.2]);

figExport(5,4,'sensitivity/overview-U')

%% function space

% plot effect of multiple params on same KPI
function multiPlot(s,fields,i_kpi,KPI)

% mapping for greek letters on y axis
yGreek = containers.Map(KPI,{'\alpha' '\beta' '\gamma' '\delta'});

k = KPI{i_kpi};

figure

for i = 1:numel(fields)
    field = fields{i};
    subplot(1,numel(fields),i)
    
    % plot effect on maximum value
    yyaxis left
    plot(s.(field).data.f,s.(field).data.([k '_max_sc']))
    ylabel(yGreek(k))
    xlabel(sprintf('(%1$s)/(%1$s)_0',field))
    
    % plot effect on z position of max
    yyaxis right
    plot(s.(field).data.f,s.(field).data.([k '_max_z'  ]))
    ylabel('z')
end

end
