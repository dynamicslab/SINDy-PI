%% This file will generate the guess of the iSINDYc left hand side.
% Last Updated: 2019/04/21
% Coded By: K
function [Data,Sym_Struct]=GuessLib(X,dX,iter,u,Highest_Poly_Order,Highest_Trig_Order,Highest_U_Order)
%% First get the size of the X matrix, determin the data length and the number of variables we have.
[Data_Length,Variable_Number]=size(X);
[~,Variable_Number_dX]=size(dX);
[~,Variable_Number_u]=size(u);
% Also create the symbolic variable
Symbol=sym('z',[Variable_Number,1]);
Symbol_dX=sym('dz',[Variable_Number,1]);
Symbol_u=sym('u',[Variable_Number_u,1]);

%% Now according the Highest Polynomial Order entered, we will calculate the data matrix.
Data=[];
Index=1;

%% First calculate the polynomial term

% Add dx term
Data(:,Index)=dX(:,4);
Sym_Struct{1,Index}=Symbol_dX(4);
Index=Index+1;


