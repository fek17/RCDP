%% CE2-03-2 Group 6
%  sensitivity

% import true constant values
constants

% set flag (so reactor.m doesn't override changed constants later)
c.Canary = true;

% sensitivity range
t.range = [0.5:0.1:1.5]; % 50 to 150%

%% start fiddling
tic
% struct c
for field_ = fieldnames(c)' 
    % convert fieldname cell to string
    field = string(field_);
    
    % check data type
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
                % convert fieldname cell to string
                k = char(k_);
                
                % save KPI data to table
                s.(field).data.([k '_max_val'])(i_fac) = o.(k).max.val;
                s.(field).data.([k '_max_z'  ])(i_fac) = o.(k).max.z;
            end
            
            % plot correlations
        end
        
        % restore original value
        c.(field) = s.(field).real;
        
    end
end
clear field_ field i_fac k_ k
toc
%% make some pretty correlation plots

t.export = true;

% for each variable
for field_ = fieldnames(s)' 
    % convert fieldname cell to string
    field = char(field_);
    
    % initialise plot
    figure
    
    % run through KPIs
    for i = 1:numel(KPI)
        % current KPI string
        k = KPI{i};
        
        % initialise subplots
        subplot(1,numel(KPI),i)
        title(['max ' KPI_latex{i}])
        xlabel(sprintf('(%1$s)/(%1$s)_0',field))

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
    figExport(18,6,['sensitivity/' field])
end
clear field_ field k_ k