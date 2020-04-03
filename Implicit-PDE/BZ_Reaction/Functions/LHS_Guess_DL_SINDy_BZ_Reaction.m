% This file will generate the necessary LHS library to identify the PDE.
% Please modify this file to code the library you want.
%
% Last Updated: 2019/06/19
% Coded By: K

function [Data,Sym_Struct]=LHS_Guess_DL_SINDy_BZ_Reaction(Data_Train,Vars,WhichEqn)
%% First get the size of the u vector.
[Data_Length,~]=size(Data_Train);

% Depending on which equation we want to solve, create different LHS guess
if WhichEqn==1
    Org_Data=Data_Train(:,1:6);
    Org_Sym=Vars(:,1:6);
elseif WhichEqn==2
    Org_Data=Data_Train(:,7:12);
    Org_Sym=Vars(:,7:12);
elseif WhichEqn==3
    Org_Data=Data_Train(:,13:18);
    Org_Sym=Vars(:,13:18);
elseif WhichEqn==4
    Org_Data=Data_Train(:,19:24);
    Org_Sym=Vars(:,19:24);
end

% Add the first term
Data=[];
Index=1;

Data(:,Index)=Org_Data(:,2);
Sym_Struct{1,Index}=Org_Sym(1,2);
Index=Index+1;

Data(:,Index)=Org_Data(:,1).*Org_Data(:,2);
Sym_Struct{1,Index}=Org_Sym(1,1)*Org_Sym(1,2);
Index=Index+1;

Data(:,Index)=Org_Data(:,1).^2.*Org_Data(:,2);
Sym_Struct{1,Index}=Org_Sym(1,1)^2*Org_Sym(1,2);
Index=Index+1;





