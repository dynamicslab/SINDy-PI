%% This function will display the ODE function.

% Last Updated: 2019/04/21
% Coded By: K

function ODEs=Print_ODEs(ODE,n_state,n_control,disp,actual)

z_vars=sym('z',[n_state;1]);
u_vars=sym('u',[n_control;1]);
d_vars=sym('dz',[n_state;1]);

if isa(ODE,'function_handle')
    if n_control~=0
        ODEs=ODE(0,z_vars,u_vars);
    else
        ODEs=ODE(0,z_vars);
    end
end

if disp==1
    if actual==1
        fprintf('\v The actual ODE of the system is/are :\n')
    else
        fprintf('\v The discovered ODE of the system is/are :\n')
    end
     for i=1:n_state
          digits(4)
          fprintf(strcat('\t',char(d_vars(i,1)),'=',char(vpa(simplify(ODEs(i,1)))),'\n'));
     end
end