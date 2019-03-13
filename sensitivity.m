%% CE2-03-2 Group 6
%  sensitivity

% import true constant values
global c
constants

% set flag (so reactor.m doesn't override changed constants later)
c.Canary = true;

% sensitivity range
t.range = [0.5:0.05:1.5]; % 50 to 150%

% set up results struct
global s
s = struct;

%% fiddle with struct c

% setup progress bar
progBar = waitbar(0,'Please wait...');

tic
fields = fieldnames(c);
n_fields = numel(fields);
for i_field = 1:n_fields
    field = string(fields{i_field});
    
    % update progress bar
    waitbar((i_field-1)/n_fields,progBar,sprintf('Evaluating field %s (%u of %u)',field,i_field,n_fields));
    
    % fiddle with number fields
    if isfloat(c.(field))
        
        doSensitivity(field,t,KPI)
    
    % fiddle with table fields
    elseif istable(c.(field))
        
        % make cell array of row and variable names character vectors
        t.rows = c.(field).Properties.RowNames;
        t.vars = c.(field).Properties.VariableNames;
        
        % replace with cell array of numbers if no names
        dim = {'rows' 'vars'};
        for i_dim = 1:2
            if numel(t.(dim{i_dim})) == 0
                t.(dim{i_dim}) = num2cell(1:size(c.(field),i_dim));
            end
        end
        
        % each variable
        for i_var = 1:numel(t.vars)
            t.var = t.vars{i_var};
            
            % each row
            for i_row = 1:numel(t.rows)
                t.row = t.rows{i_row};
                if isfloat(t.row)
                    t.row_ = num2str(t.row);
                else
                    t.row_ = t.row;
                end
                
                % label for this field
                t.table_field = [t.var '_' t.row_];
                
                doSensitivity(field,t,KPI)
            
            end
        end
        
    end
    
end
clear i_field field_ field i_fac k_ k dim i_dim i_var i_row
toc

% close progress bar
close(progBar)

%% make some pretty correlation plots

% export?
global printFlag
printFlag = false;

% for each variable
%for field_ = fieldnames(s)'
for field_ = {'P0' 'T0' 'U'}
    field = char(field_);
    
    % initialise plot
    figure
    
    % run through KPIs
    for i_kpi = 1:numel(KPI)
        % current KPI string
        k = KPI{i_kpi};
        
        % initialise subplots
        subplot(2,ceil(numel(KPI)/2),i_kpi)
        title(KPI_latex{i_kpi})
        xlabel(sprintf('(%1$s)/(%1$s)_0',strrep(field,'_','\_')))

        % calculate scaled maximum value
        s.(field).data.([k '_max_sc']) = s.(field).data.([k '_max_val']) ./ dot(s.(field).data.([k '_max_val']),s.(field).data.f==1);
        
        % plot effect on maximum value
        yyaxis left
        plot(s.(field).data.f,s.(field).data.([k '_max_sc']))
        ylabel('(max)/(max)_0')
        
        % plot effect on z position of max
        yyaxis right
        plot(s.(field).data.f,s.(field).data.([k '_max_z'  ]))
        ylabel('z_{max} / m')
    end
    
    % export
    figExport(12,12,['sensitivity/' field])
end
clear field_ field i_kpi k_ k

%% function space

function doSensitivity(field,t,KPI)
global c v s o
o = struct;
v = struct;

% define field label
if isfloat(c.(field))
    field_label = field;
elseif istable(c.(field))
    field_label = t.table_field;
end

% store original value
s.(field_label).real = readNestedField(field,t);

% repeat for each multiplication factor
for i_fac = 1:numel(t.range)
    
    % initialise results table on first run
    if i_fac == 1
        % empty columns
        s.(field_label).data = array2table(zeros(numel(t.range),2*numel(KPI)+1));
        s.(field_label).data.Properties.VariableNames = [{'f'} strcat(KPI,'_max_val') strcat(KPI, '_max_z')];
        
        % save multiplication factors
    	s.(field_label).data.f = t.range';
    end
    
    % get current multiplication factor
	t.r = t.range(i_fac);
    
    % set new value of constant
    newVal = s.(field_label).real * t.r;
    writeNestedField(field,t,newVal)
    
    % calculate new KPIs
    try
        % fun with scoping variables
        i = 0;
        j = 0;
        n_j = "";
        
        % run the reactor code
        reactor
        
    catch ME % catch errors
        s.(field_label).error = ME;
        getReport(ME);
    end
    
    % save results for each KPI
    for k_ = fieldnames(o)'
        k = char(k_);
        
        % save KPI data to table
        s.(field_label).data.([k '_max_val'])(i_fac) = o.(k).max.val;
        s.(field_label).data.([k '_max_z'  ])(i_fac) = o.(k).max.z;
    end
    
end

% restore original value
origVal = s.(field_label).real;
writeNestedField(field,t,origVal)

% read the field
function out = readNestedField(field,t)
if isfloat(c.(field))
    out = c.(field);
elseif istable(c.(field))
    out = c.(field).(t.var)(t.row);
end
end

% write the field
function writeNestedField(field,t,in)
if isfloat(c.(field))
    c.(field) = in;
elseif istable(c.(field))
    c.(field).(t.var)(t.row) = in;
end
end

end
