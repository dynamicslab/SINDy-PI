%% This function will calculate the prediction error between the discovered system and actual test data
% There are two methods to calculate the prediction error.
%
% Last Updated: 2019/06/20
% Coded By: K
%%
function [Score]=Get_Score(LHS_Data_Test,Data_Test,name)

% First create a function handel for the PDE file
PDE_func=str2func(strcat('@(z)',name,'(z)'));

% Then get the LHS estimation
LHS_Es=PDE_func(Data_Test);

% Now calculate the l2 norm
Score=norm(LHS_Data_Test-LHS_Es)/norm(LHS_Data_Test);

end



