% 2019 Dummy Test
% File name PA_Reactor_2019 MATLAB 2018b
% 10s Feed Gas inlet temperature goes up 1K.
% 20s Feed Gas inlet temperature goes up 1K
% 30s Feed Gas inlet temperature goes down 2K
% 50s Feed Gas inlet temperature goes down 2K, Cooling medium inlet temperature goes up 3K
% 60s simulations end.

load_system('PA_Reactor_2019');

set_param('PA_Reactor_2019/Reactants Inlet temperature  (K)', 'OutValues', "[690 691 692 692 692 690].'", 'tsamp', '10');
set_param('PA_Reactor_2019/Cooling fluid inlet temperature (K)', 'OutValues', "[680 680 678 677 680 680].'", 'tsamp', '10');

set_param('PA_Reactor_2019', 'Solver', 'ode15s', 'StopTime', num2str(60));

sim('PA_Reactor_2019')