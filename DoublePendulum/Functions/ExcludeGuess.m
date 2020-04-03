%% This function will exclude the left hand side guess from the right hand side library.
% Last Updated: 2019/04/22
% Coded By: K
function [RHS_Data,RHS_Struct]=ExcludeGuess(SINDy_Data,SINDy_Struct,LHS_Sym)
% First get the size of the library
[n,m]=size(SINDy_Data);
% Then run a for loop and generate a new library
j=1;
for i=1:m
    if SINDy_Struct{i}~=LHS_Sym
        RHS_Data(:,j)=SINDy_Data(:,i);
        RHS_Struct(1,j)=SINDy_Struct(1,i);
        j=j+1;
    else 
        % For use of debuging 
        SINDy_Struct{i};
    end
end















