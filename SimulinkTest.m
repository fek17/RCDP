%% Tests for Simulink Temperature Control Loop

%% Default test
set_param('PA_Reactor_2019/Reactants Inlet temperature  (K)', 'OutValues', "[690 691 692 692 692 690]", 'tsamp', '10');
set_param('PA_Reactor_2019/Cooling fluid inlet temperature (K)', 'OutValues', "[680 680 678 677 680 680]", 'tsamp', '10');

%% Alessandro test
set_param('PA_Reactor_2019/Reactants Inlet temperature  (K)', 'OutValues', "[]", 'tsamp', '10');
set_param('PA_Reactor_2019/Cooling fluid inlet temperature (K)', 'OutValues', "[680 680 685 690 695 705 705 705 705 705]", 'tsamp', '30');

%% Random test

% reactant
T(1).lims = [665 715];
T(1).step = 10;
T(1).block = 'PA_Reactor_2019/Reactants Inlet temperature  (K)';

% coolant
T(2).lims = [665 695];
T(2).step = 30;
T(2).block = 'PA_Reactor_2019/Cooling fluid inlet temperature (K)';

for i = 1:numel(T)
    
    % map expected limits to 4*sigma (~4% chance of falling outside)
    T(i).sigma = ( T(i).lims(2) - T(i).lims(1) )/4;
    
    % map midpoint of limits to mu
    T(i).mu = ( T(i).lims(1) + T(i).lims(2) )/2;
    
    % make array of random integers (normal distribution)
    T(i).array = randn(1,10) * T(i).sigma + T(i).mu;
    
    % format array into string
    T(i).string = join(['[' compose("%u",T(i).array) ']']);
    
    % apply parameters
    set_param(T(i).block, 'OutValues', T(i).string, 'tsamp', num2str(T(i).step));

end