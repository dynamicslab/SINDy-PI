%% This function will generate the PDE function based on calculation result of the symbolic expression.
% This file has been modified to perform the parallel computing.
% Last Update: 2019/06/20
% Coded By: K

function Generate_PDE_RHS(Sindy_PDEs,Vars,name)
f_x=Sindy_PDEs;
matlabFunction(f_x,'File',name,'Optimize',true,'Vars',{Vars});










