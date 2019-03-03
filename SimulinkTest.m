% Tests for Simulink Temperature Control Loop
load_system('PA_Reactor_2019');
set_param('PA_Reactor_2019/Cooling fluid inlet temperature (K)', 'OutValues', "[680 680 685 690 695 705 705 705 705 705].'", 'tsamp', '30');
set_param('PA_Reactor_2019', 'Solver', 'ode15s', 'StopTime', num2str(60));
sim('PA_Reactor_2019')