% This file will generate the necessary RHS library to identify the PDE.
% Please modify this file to code any library you want.
%
% Last Updated: 2019/06/19
% Coded By: K

function [Data,Sym_Struct]=SINDyLib_DL_SINDy_BZ_Reaction(Data_Train,Vars,WhichEqn)
%% First get the size of the vector.
[Data_Length,~]=size(Data_Train);

Data=[];
Index=1;

% Add constant terms in your library
Data(:,1)=ones(Data_Length,1);
Sym_Struct{1,1}=1;
Index=Index+1;

% Depending on which equation we want to solve, create different LHS guess
if WhichEqn==1
    Org_Data=Data_Train(:,1:6);
    Org_Sym=Vars(:,1:6);
    Sup_Data(:,1:3)=[Data_Train(:,7) Data_Train(:,13) Data_Train(:,19)];
    Sup_Sym(1,1:3)=[Vars(:,7) Vars(:,13) Vars(:,19)];
elseif WhichEqn==2
    Org_Data=Data_Train(:,7:12);
    Org_Sym=Vars(:,7:12);   
    Sup_Data(:,1:3)=[Data_Train(:,1) Data_Train(:,13) Data_Train(:,19)];
    Sup_Sym(1,1:3)=[Vars(:,1) Vars(:,13) Vars(:,19)];
elseif WhichEqn==3
    Org_Data=Data_Train(:,13:18);
    Org_Sym=Vars(:,13:18);   
    Sup_Data(:,1:3)=[Data_Train(:,1) Data_Train(:,7) Data_Train(:,19)];
    Sup_Sym(1,1:3)=[Vars(:,1) Vars(:,7) Vars(:,19)];
elseif WhichEqn==4
    Org_Data=Data_Train(:,19:24);
    Org_Sym=Vars(:,19:24);
    Sup_Data(:,1:3)=[Data_Train(:,1) Data_Train(:,7) Data_Train(:,13)];
    Sup_Sym(1,1:3)=[Vars(:,1) Vars(:,7) Vars(:,13)];
end

for pin=1:3
    Data(:,Index)=Org_Data(:,1).^pin;
    Sym_Struct{1,Index}=Org_Sym(1,1).^pin;
    Index=Index+1;
end

for pin=1:3
    Data(:,Index)=Sup_Data(:,pin);
    Sym_Struct{1,Index}=Sup_Sym(1,pin);
    Index=Index+1;
end

for pin=1:3
    Data(:,Index)=Org_Data(:,1).*Sup_Data(:,pin);
    Sym_Struct{1,Index}=Org_Sym(1,1)*Sup_Sym(1,pin);
    Index=Index+1;
end

Data(:,Index)=Org_Data(:,2);
Sym_Struct{1,Index}=Org_Sym(1,2);
Index=Index+1;

Data(:,Index)=Org_Data(:,1).*Org_Data(:,2);
Sym_Struct{1,Index}=Org_Sym(1,1)*Org_Sym(1,2);
Index=Index+1;

for pin=1:4
    Data(:,Index)=Org_Data(:,2+pin);
    Sym_Struct{1,Index}=Org_Sym(1,2+pin);
    Index=Index+1;
end

for pin=1:4
    Data(:,Index)=Org_Data(:,1).*Org_Data(:,2+pin);
    Sym_Struct{1,Index}=Org_Sym(:,1)*Org_Sym(1,2+pin);
    Index=Index+1;
end











