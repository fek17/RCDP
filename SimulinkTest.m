%% Tests for Simulink Temperature Control Loop

%% Default test
set_param('PA_Reactor_2019/Reactants Inlet temperature  (K)', 'OutValues', "[690 691 692 692 692 690]", 'tsamp', '10');
set_param('PA_Reactor_2019/Cooling fluid inlet temperature (K)', 'OutValues', "[680 680 678 677 680 680]", 'tsamp', '10');

%% Alessandro test
set_param('PA_Reactor_2019/Reactants Inlet temperature  (K)', 'OutValues', "[]", 'tsamp', '10');
set_param('PA_Reactor_2019/Cooling fluid inlet temperature (K)', 'OutValues', "[680 680 685 690 695 705 705 705 705 705]", 'tsamp', '30');
