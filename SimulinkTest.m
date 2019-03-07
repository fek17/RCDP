%% Tests for Simulink Temperature Control Loop

%% Default test
set_param('PA_Reactor_2019/Reactants Inlet temperature  (K)', 'OutValues', "[690 691 692 692 692 690]", 'tsamp', '10');
set_param('PA_Reactor_2019/Cooling fluid inlet temperature (K)', 'OutValues', "[680 680 678 677 680 680]", 'tsamp', '10');

%% Alessandro test
set_param('PA_Reactor_2019/Reactants Inlet temperature  (K)', 'OutValues', "[]", 'tsamp', '10');
set_param('PA_Reactor_2019/Cooling fluid inlet temperature (K)', 'OutValues', "[680 680 685 690 695 705 705 705 705 705]", 'tsamp', '30');

%% Random test

% reactant
T(1).range = [660 720];
T(1).step = 5;
T(1).block = 'PA_Reactor_2019/Reactants Inlet temperature  (K)';

% coolant
T(2).range = [660 700];
T(2).step = 10;
T(2).block = 'PA_Reactor_2019/Cooling fluid inlet temperature (K)';

for i = 1:numel(T)
    
    % make array of random integers
    T(i).array = randi(T(i).range,1,10);
    
    % format array into string
    T(i).string = join(['[' compose("%u",T(i).array) ']']);
    
    % apply parameters
    set_param(T(i).block, 'OutValues', T(i).string, 'tsamp', num2str(T(i).step));

end