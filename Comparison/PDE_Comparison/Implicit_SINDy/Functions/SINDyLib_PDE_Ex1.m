% This file will generate the necessary RHS library to identify the PDE.
% Please modify this file to code the library you want.
%
% Last Updated: 2019/06/19
% Coded By: K

function [Data,Sym_Struct]=SINDyLib_PDE_Ex1(U,Ut,Utt,Ux,Uxx,Uxxx)
%% First get the size of the u vector.
[Data_Length,~]=size(U);

%Also create the symbolic variable
syms u ut utt ux uxx uxxx

Data=[];
Index=1;

Org_Data=[U Ux Uxx Uxxx];
Org_Sym=[u ux uxx uxxx];

% Add terms in your library
Data(:,1)=ones(Data_Length,1);
Sym_Struct{1,1}=1;
Index=Index+1;

Data(:,Index)=Ut;
Sym_Struct{1,Index}=ut;
Index=Index+1;

Data(:,Index)=U.*Ut;
Sym_Struct{1,Index}=u*ut;

for i=1:size(Org_Data,2)
    Index=Index+1;
    Data(:,Index)=Org_Data(:,i);
    Sym_Struct{1,Index}=Org_Sym(1,i);
end

for i=2:size(Org_Data,2)
    Index=Index+1;
    Data(:,Index)=Org_Data(:,1).*Org_Data(:,i);
    Sym_Struct{1,Index}=Org_Sym(1,1)*Org_Sym(1,i);
end





