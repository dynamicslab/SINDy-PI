% This library is coded for the double pendulum example
% Last Updated: 2019/07/30
% Coded By: K

function [Data,Sym_Struct]=SINDyLib(X,dX,iter,u,Highest_Poly_Order,Highest_Trig_Order,Highest_U_Order,Highest_dPoly_Order)
%% First get the size of the X matrix, determin the data length and the number of variables we have.
[Data_Length,Variable_Number]=size(X);
[~,Variable_Number_dX]=size(dX);
[~,Variable_Number_u]=size(u);
%Also create the symbolic variable
Symbol=sym('z',[Variable_Number,1]);
Symbol_dX=sym('dz',[Variable_Number,1]);
Symbol_u=sym('u',[Variable_Number_u,1]);

%Now according the Highest Polynomial Order entered, we will calculate the data matrix.
Data=[];
Index=1;

%% First calculate the polynomial term
%Order zero:
% Data(:,Index)=ones(Data_Length,1);
% Sym_Struct{1,Index}=1;

% Form basis vector
Basis=[X(:,1) X(:,2) X(:,1)-X(:,2) X(:,1)+X(:,2)];
Basis_Sym=[Symbol(1) Symbol(2) Symbol(1)-Symbol(2) Symbol(1)+Symbol(2)];

% Add the trigonometric form
for i=1:size(Basis,2)
   Data(:,Index)=sin(Basis(:,i));
   Sym_Struct{1,Index}=sin(Basis_Sym(1,i));
   Index=Index+1;
end

for i=1:size(Basis,2)
   Data(:,Index)=cos(Basis(:,i));
   Sym_Struct{1,Index}=cos(Basis_Sym(1,i));
   Index=Index+1;
end

% Add dx term
pin1=Index;
for i=1:4
    Data(:,Index)=dX(:,i);
    Sym_Struct{1,Index}=Symbol_dX(i);
    Index=Index+1;
end

for i=1:4
    Data(:,Index)=dX(:,i).^2;
    Sym_Struct{1,Index}=Symbol_dX(i)^2;
    Index=Index+1;
end
pin2=Index;

% Add the combanition of two
for i=pin1:pin2-1
    for j=1:pin1-1
        Data(:,Index)=Data(:,i).*Data(:,j);
        Sym_Struct{1,Index}=Sym_Struct{1,i}*Sym_Struct{1,j};
        Index=Index+1;
    end
end

