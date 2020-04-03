% This library is coded for the single pendulum example
% Last Updated: 2019/04/22
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
Index=1;
Data(:,Index)=ones(Data_Length,1);
Sym_Struct{1,Index}=1;

%Order One:
if Highest_Poly_Order>=1
    for i=1:Variable_Number
        Index=Index+1;
        Data(:,Index)=X(:,i);
        Sym_Struct{1,Index}=Symbol(i,1);
    end
end
pin1=Index;

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
                for pi=1:Variable_Number
                    Index=Index+1;
                    Data(:,Index)=X(:,i).*X(:,j).*X(:,k).*X(:,pi);
                    Sym_Struct{1,Index}=Symbol(i,1)*Symbol(j,1)*Symbol(k,1)*Symbol(pi,1);
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
    %  
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
pin1=Index;

%Order Three:
if Highest_Trig_Order>=3
    for i=1:Variable_Number
        for j=1:Variable_Number
            Index=Index+1;
            Data(:,Index)=X(:,i).^2.*sin(X(:,j)).*cos(X(:,j));
            Sym_Struct{1,Index}=Symbol(i,1)^2*sin(Symbol(j,1))*cos(Symbol(j,1));
        end
    end
end

%Order Four:
if Highest_Trig_Order>=4
    for i=1:Variable_Number
        for j=1:Variable_Number
            Index=Index+1;
            Data(:,Index)=X(:,i).^2.*sin(X(:,j));
            Sym_Struct{1,Index}=Symbol(i,1)^2*sin(Symbol(j,1));
        end
    end
    %
    for i=1:Variable_Number
        for j=1:Variable_Number
            Index=Index+1;
            Data(:,Index)=X(:,i).^2.*cos(X(:,j));
            Sym_Struct{1,Index}=Symbol(i,1)^2*cos(Symbol(j,1));
        end
    end
    %
    for i=1:Variable_Number
        for j=1:Variable_Number
            Index=Index+1;
            Data(:,Index)=sin(X(:,i)).*cos(X(:,j));
            Sym_Struct{1,Index}=sin(Symbol(i,1))*cos(Symbol(j,1));
        end
    end
end

pin=Index;

%% From here, we add the dX*Theta elements in our data.
if Highest_dPoly_Order>=1
    for k=1:pin1
        Index=Index+1;
        Data(:,Index)=dX(:,1).*Data(:,k);
        Sym_Struct{1,Index}=Symbol_dX(iter,1)*(Sym_Struct{1,k});
    end
end

%% Frome here, we add the u*Theta elements in our data
if Highest_U_Order>=1
    for j=1:Variable_Number_u
        for k=1:pin
            Index=Index+1;
            Data(:,Index)=u(:,j).*Data(:,k);
            Sym_Struct{1,Index}=Symbol_u(j,1)*(Sym_Struct{1,k});
        end
    end
else
    for j=1:Variable_Number_u
        Index=Index+1;
        Data(:,Index)=u(:,j);
        Sym_Struct{1,Index}=Symbol_u(j,1);
    end
end

