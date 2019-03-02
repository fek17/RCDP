%% CE2-03-2 Group 6
%  sensitivity

% import true constant values
constants

% sensitivity range
t.S_range = [0.5:0.5:1.5]; % 50 to 150% 

fiddlerMain(t.S_range)
fiddlerExpander({'f' 'n0'})
function fiddlerMain(r)

% inputs
%   r   (vector) range of multiplication factors for test

global c

% iterate through struct c
for thisField = fieldnames(c)'
    % pass through to sorter (since structs can contain anything)
    %fiddlerSort(thisField{1},r)
end
end

function fiddlerSort(x,r)

% pass x on to a type-specific fiddler function
if class(c.(x)) == 'double'
    fiddlerDouble(x,r)
elseif class(c.(x)) == 'struct'
    fiddlerStruct(x,r)
elseif class(c.(x)) == 'table'
    fiddlerTable(x,r)
end

end

function fiddlerDouble(x,r)
% run through for each multiplication factor
for f = r
    % calculate new KPIs
    % plot correlations
end
end

function fiddlerStruct(x,r)
% iterate across all fields of struct
for thisField = fieldnames(c)
    % pass through to sorter (since structs can contain anything)
    fiddlerSort(x.(cellstr(thisField)),r)
end
end

function fiddlerTable(x,r)
% iterate across cells of table
end

function fiddlerExpander(x)

% inputs
%   x   (cell) layered struct field names, excl. base c
global c
if numel(x) == 1
    c.(x{1})
elseif numel(x) == 2
    c.(x{1}).(x{2})
elseif numel(x) == 3
    c.(x{1}).(x{2}).(x{3})
end

end