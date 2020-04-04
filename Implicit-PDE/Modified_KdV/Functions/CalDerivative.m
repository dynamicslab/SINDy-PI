%% This file will generate the first derivative of the data.
% Coded By: K
% Last Updated: 2019/06/19
%% Function file
function Dev=CalDerivative(x,dx,d)
% First we get the information of the data length. The x should be a nx1
% vector.
[n,m]=size(x);

[len,Index]=max([n,m]);

if Index==1
    Dev=zeros(len,1);
else
    Dev=zeros(1,len);
end

%% Define the coeficient for different orders of derivative
if d==1
    p1=1/12;p2=-2/3;p3=0;p4=2/3;p5=-1/12;
elseif d==2
    p1=-1/12;p2=4/3;p3=-5/2;p4=4/3;p5=-1/12;
elseif d==3
    p1=-1/2;p2=1;p3=0;p4=-1;p5=1/2;
end

% Calculate the derivative of the middel point
for i=3:len-2
    Dev(i)=(p1*x(i-2)+p2*x(i-1)+p3*x(i)+p4*x(i+1)+p5*x(i+2));
    if d==1
        Dev(i)=Dev(i)/dx;
    elseif d==2
        Dev(i)=Dev(i)/dx^2;
    elseif d==3
        Dev(i)=Dev(i)/dx^3;
    end
end

%% Ge the derivative of first two points using forward difference
if d==1
    q1=-3/2;q2=2;q3=-1/2;q4=0;q5=0;
elseif d==2
    q1=2;q2=-5;q3=4;q4=-1;q5=0;
elseif d==3
    q1=-5/2;q2=9;q3=-12;q4=7;q5=-3/2;
end

for i=1:2
    Dev(i)=(q1*x(i)+q2*x(i+1)+q3*x(i+2)+q4*x(i+3)+q5*x(i+4));
    if d==1
        Dev(i)=Dev(i)/dx;
    elseif d==2
        Dev(i)=Dev(i)/dx^2;
    elseif d==3
        Dev(i)=Dev(i)/dx^3;
    end
end

%% Ge the derivative of last two points using backward difference
if d==1
    m1=3/2;m2=-2;m3=1/2;m4=0;m5=0;
elseif d==2
    m1=2;m2=-5;m3=4;m4=-1;m5=0;
elseif d==3
    m1=5/2;m2=-9;m3=12;m4=-7;m5=3/2;
end

for i=len-1:len
    Dev(i)=(m1*x(i)+m2*x(i-1)+m3*x(i-2)+m4*x(i-3)+m5*x(i-4));
    if d==1
        Dev(i)=Dev(i)/dx;
    elseif d==2
        Dev(i)=Dev(i)/dx^2;
    elseif d==3
        Dev(i)=Dev(i)/dx^3;
    end
end







