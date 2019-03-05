%% CE2-03-2 Group 6
%  sensitivity

% import true constant values
constants

% set flag (so reactor.m doesn't override changed constants later)
c.Canary = true;

% sensitivity range
t.range = [0.5:0.1:1.5]; % 50 to 150%

% set up structs
s = struct;

%% fiddle with struct c
tic
for field_ = fieldnames(c)' 
    field = string(field_);
    
    % fiddle with number fields
    if isfloat(c.(field))
        
        % store original value
        s.(field).real = c.(field);
        
        % repeat for each multiplication factor
        for i_fac = 1:numel(t.range)
            
            % initialise results table on first run
            if i_fac == 1
                s.(field).data = array2table(zeros(numel(t.range),2*numel(KPI)+1));
                s.(field).data.Properties.VariableNames = [{'f'} strcat(KPI,'_max_val') strcat(KPI, '_max_z')];
            end
            
            % get current multiplication factor & save
            t.r = t.range(i_fac);
            s.(field).data{i_fac,'f'} = t.r;
            
            % set new value of constant
            c.(field) = s.(field).real * t.r;
            
            % calculate new KPIs
            try
                reactor
            catch ME % catch errors
                s.(field).error = ME;
            end
            
            % save results for each KPI
            for k_ = fieldnames(o)'
                k = char(k_);
                
                % save KPI data to table
                s.(field).data.([k '_max_val'])(i_fac) = o.(k).max.val;
                s.(field).data.([k '_max_z'  ])(i_fac) = o.(k).max.z;
            end
            
        end
        
        % restore original value
        c.(field) = s.(field).real;
    
    end
    
    % fiddle with table fields
    if istable(c.(field))
        
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
            var = t.vars{i_var};
            
            % each row
            for i_row = 1:numel(t.rows)
                row = t.rows{i_row};
                if isfloat(row)
                    row_ = num2str(row);
                else
                    row_ = row;
                end
                
                % store original value
                s.([var '_' row_]).real = c.(field){row,var};
                
                % repeat for each multiplication factor
                for i_fac = 1:numel(t.range)

                    % initialise results table on first run
                    if i_fac == 1
                        s.([var '_' row_]).data = array2table(zeros(numel(t.range),2*numel(KPI)+1));
                        s.([var '_' row_]).data.Properties.VariableNames = [{'f'} strcat(KPI,'_max_val') strcat(KPI, '_max_z')];
                    end

                    % get current multiplication factor & save
                    t.r = t.range(i_fac);
                    s.([var '_' row_]).data{i_fac,'f'} = t.r;

                    % set new value of constant
                    c.(field){row,var} = s.([var '_' row_]).real * t.r;

                    % calculate new KPIs
                    try
                        reactor
                    catch ME % catch errors
                        s.([var '_' row_]).error = ME;
                    end

                    % save results for each KPI
                    for k_ = fieldnames(o)'
                        % convert fieldname cell to string
                        k = char(k_);

                        % save KPI data to table
                        s.([var '_' row_]).data.([k '_max_val'])(i_fac) = o.(k).max.val;
                        s.([var '_' row_]).data.([k '_max_z'  ])(i_fac) = o.(k).max.z;
                    end

                end

                % restore original value
                c.(field){row,var} = s.([var '_' row_]).real;
            
            end
        end
        
    end
    
end
clear field_ field i_fac k_ k dim i_dim i_var var i_row row row_
toc

%% make some pretty correlation plots

t.export = true;

% for each variable
for field_ = fieldnames(s)' 
    field = char(field_);
    
    % initialise plot
    figure
    
    % run through KPIs
    for i = 1:numel(KPI)
        % current KPI string
        k = KPI{i};
        
        % initialise subplots
        subplot(2,ceil(numel(KPI)/2),i)
        title(['max ' KPI_latex{i}])
        xlabel(sprintf('(%1$s)/(%1$s)_0',strrep(field,'_','\_')))

        % calculate scaled maximum value
        s.(field).data.([k '_max_sc']) = s.(field).data.([k '_max_val']) ./ dot(s.(field).data.([k '_max_val']),s.(field).data.f==1);
        
        % plot effect on maximum value
        yyaxis left
        plot(s.(field).data.f,s.(field).data.([k '_max_sc']))
        ylabel('(value)/(value)_0')
        
        % plot effect on z position of max
        yyaxis right
        plot(s.(field).data.f,s.(field).data.([k '_max_z'  ]))
        ylabel('z / m')
    end
    
    % export
    figExport(12,12,['sensitivity/' field])
end
clear field_ field k_ k
