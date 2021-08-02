% This library is coded for the double pendulum example
% Last Updated: 2019/07/30
% Coded By: K.Kahirman

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

% Form basis vector
Basis=[X(:,1) X(:,2) X(:,1)-X(:,2) X(:,1)-2*X(:,2) 2*X(:,1)-X(:,2) 2*X(:,1)-2*X(:,2)];
Basis_Sym=[Symbol(1) Symbol(2) Symbol(1)-Symbol(2) Symbol(1)-2*Symbol(2) 2*Symbol(1)-Symbol(2) 2*Symbol(1)-2*Symbol(2)];

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

for i=3:size(Basis,2)
   Data(:,Index)=cos(Basis(:,i)).^2;
   Sym_Struct{1,Index}=cos(Basis_Sym(1,i))^2;
   Index=Index+1;
end

% Adding following terms will reduce the noise robustness
% for i=3:size(Basis,2)
%    Data(:,Index)=sin(Basis(:,1)).*cos(Basis(:,i));
%    Sym_Struct{1,Index}=sin(Basis_Sym(1,1))*cos(Basis_Sym(1,i));
%    Index=Index+1;
% end
% 
% for i=3:size(Basis,2)
%    Data(:,Index)=cos(Basis(:,1)).*cos(Basis(:,i));
%    Sym_Struct{1,Index}=cos(Basis_Sym(1,1))*cos(Basis_Sym(1,i));
%    Index=Index+1;
% end
% 
% for i=3:size(Basis,2)
%    Data(:,Index)=sin(Basis(:,2)).*cos(Basis(:,i));
%    Sym_Struct{1,Index}=sin(Basis_Sym(1,2))*cos(Basis_Sym(1,i));
%    Index=Index+1;
% end
% 
% for i=3:size(Basis,2)
%    Data(:,Index)=cos(Basis(:,2)).*cos(Basis(:,i));
%    Sym_Struct{1,Index}=cos(Basis_Sym(1,2))*cos(Basis_Sym(1,i));
%    Index=Index+1;
% end

% Add polynomial term
Data(:,Index)=X(:,3);
Sym_Struct{1,Index}=Symbol(3);
Index=Index+1;

Data(:,Index)=X(:,4);
Sym_Struct{1,Index}=Symbol(4);
Index=Index+1;

for i=3:size(Basis,2)
   Data(:,Index)=X(:,3).*sin(Basis(:,i));
   Sym_Struct{1,Index}=Symbol(3)*sin(Basis_Sym(1,i));
   Index=Index+1;
end

for i=3:size(Basis,2)
   Data(:,Index)=X(:,4).*sin(Basis(:,i));
   Sym_Struct{1,Index}=Symbol(4)*sin(Basis_Sym(1,i));
   Index=Index+1;
end

for i=3:size(Basis,2)
   Data(:,Index)=X(:,3).^2.*sin(Basis(:,i));
   Sym_Struct{1,Index}=Symbol(3)^2*sin(Basis_Sym(1,i));
   Index=Index+1;
end

for i=3:size(Basis,2)
   Data(:,Index)=X(:,4).^2.*sin(Basis(:,i));
   Sym_Struct{1,Index}=Symbol(4)^2*sin(Basis_Sym(1,i));
   Index=Index+1;
end

% Add dx term
Data(:,Index)=dX(:,1);
Sym_Struct{1,Index}=Symbol_dX(iter);
Index=Index+1;

for i=3:size(Basis,2)
   Data(:,Index)=dX(:,1).*cos(Basis(:,i)).^2;
   Sym_Struct{1,Index}=Symbol_dX(iter)*cos(Basis_Sym(1,i)).^2;
   Index=Index+1;
end

% Add constant
%Order zero:
Data(:,Index)=ones(Data_Length,1);
Sym_Struct{1,Index}=1;



