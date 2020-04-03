%% This file will generate the guess of the iSINDYc left hand side.
% Last Updated: 2019/04/21
% Coded By: K
function [P_Data,P_sym]=GuessLib(X,dX,iter,u,Highest_Poly_Order,Highest_Trig_Order,Highest_U_Order)
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

%Order Zero:
Data(:,Index)=ones(Data_Length,1);
Sym_Struct{1,Index}=sym(1);

%Order One:
if Highest_Poly_Order>=1
    for i=1:Variable_Number
        Index=Index+1;
        Data(:,Index)=X(:,i);
        Sym_Struct{1,Index}=Symbol(i,1);
    end
end

%Order Two:
if Highest_Poly_Order>=2
    for i=1:Variable_Number
        for j=1:Variable_Number
            Index=Index+1;
            Data(:,Index)=X(:,i).*X(:,j);
            Sym_Struct{1,Index}=Symbol(i,1)*Symbol(j,1);
        end
    end
end

%Order Three:
if Highest_Poly_Order>=3
    for i=1:Variable_Number
        for j=1:Variable_Number
            for k=1:Variable_Number
                Index=Index+1;
                Data(:,Index)=X(:,i).*X(:,j).*X(:,k);
                Sym_Struct{1,Index}=Symbol(i,1)*Symbol(j,1)*Symbol(k,1);
            end
        end
    end
end

%Order Four:
if Highest_Poly_Order>=4
    for i=1:Variable_Number
        for j=1:Variable_Number
            for k=1:Variable_Number
                for p=1:Variable_Number
                    Index=Index+1;
                    Data(:,Index)=X(:,i).*X(:,j).*X(:,k).*X(:,p);
                    Sym_Struct{1,Index}=Symbol(i,1)*Symbol(j,1)*Symbol(k,1)*Symbol(p,1);
                end
            end
        end
    end
end

%% Then add the Trigonometric Function in the output data:

%Order One:
if Highest_Trig_Order>=1
    for i=1:Variable_Number
        Index=Index+1;
        Data(:,Index)=sin(X(:,i));
        Sym_Struct{1,Index}=sin(Symbol(i,1));
    end
    for i=1:Variable_Number
        Index=Index+1;
        Data(:,Index)=cos(X(:,i));
        Sym_Struct{1,Index}=cos(Symbol(i,1));
    end
end

%Order Two:
if Highest_Trig_Order>=2
%     for i=1:Variable_Number
%         Index=Index+1;
%         Data(:,Index)=sin(X(:,i)).^2;
%         Sym_Struct{1,Index}=sin(Symbol(i,1))^2;     
%     end
    %
    for i=1:Variable_Number
        Index=Index+1;
        Data(:,Index)=cos(X(:,i)).^2;
        Sym_Struct{1,Index}=cos(Symbol(i,1))^2;     
    end
end

pin=Index;
j=0;
%% 
for k=1:Variable_Number_dX
    j=j+1;
    P_Data(:,j)=dX(:,1);
    P_sym{1,j}=Symbol_dX(iter,1);
end
%
for i=1:Variable_Number_u
    for k=1:pin
        j=j+1;
        P_Data(:,j)=u(:,i).*Data(:,k);
        P_sym{1,j}=Symbol_u(i,1)*(Sym_Struct{1,k});
    end
end


