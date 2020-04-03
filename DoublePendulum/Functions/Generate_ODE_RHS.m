%% This function will generate the ODE function based on calculation result of the sybolic expression
% Last Update: 2019/04/22
% Coded By: K

function Generate_ODE_RHS(Sindy_ODEs,var_num_state,var_num_control)
z=sym('z',[var_num_state,1]);
u=sym('u',[var_num_control,1]);
syms t
f= matlabFunction(Sindy_ODEs,'File','Sindy_ODE_RHS','Optimize',true,'Vars',{t,z,u});










