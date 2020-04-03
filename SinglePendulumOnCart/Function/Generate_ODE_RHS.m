%% This function will generate the ODE function based on calculation result of the sybolic expression
% Last Update: 2019/04/22
% Coded By: K
%(Modified so that matrix demension match)

function Generate_ODE_RHS(Sindy_ODEs,var_num_state,var_num_control)
z=sym('z',[1,var_num_state]);
u=sym('u',[1,var_num_control]);
syms t
f= matlabFunction(Sindy_ODEs,'File','Sindy_ODE_RHS','Optimize',true,'Vars',{t,z,u});










