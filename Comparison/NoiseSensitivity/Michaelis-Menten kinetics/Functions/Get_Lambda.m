%% This function is used to get the values of the lambda. This function is necessary to perform the parallel computing.
% Coded By: K
% Last Updated: 2019/05/14
%%
function lambda=Get_Lambda(j)
% Update the lambda
if j<=10
    % Reset the lambda
    dlambda=(0.1-0)/10;
    % Set the lamdba ranging from 0 ~ 0.1
   lambda=j*dlambda;
elseif j<=28
    % Reset the lambda
    dlambda=(1-0.1)/18;
    % Set the lambda ranging from 0.1 ~ 1
   lambda=(j-10)*dlambda+0.1;
elseif j<=68
    % Reset the lambda
    dlambda=(5-1)/40;
    % Set the lambda ranging from 1 ~ 5
   lambda=(j-28)*dlambda+1;
end