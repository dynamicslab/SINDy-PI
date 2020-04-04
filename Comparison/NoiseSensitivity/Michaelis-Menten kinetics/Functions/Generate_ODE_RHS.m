%% This function will generate the ODE function based on calculation result of the symbolic expression.
% This file has been modified to perform the parallel computing.
% Last Update: 2019/05/15
% Coded By: K
%%
function Generate_ODE_RHS(Sindy_ODEs,var_num_state,var_num_control,name)
z=sym('z',var_num_state);
u=sym('u',var_num_control);
syms t
f= matlabFunction(Sindy_ODEs,'File',name,'Optimize',true,'Vars',{t,z,u});










