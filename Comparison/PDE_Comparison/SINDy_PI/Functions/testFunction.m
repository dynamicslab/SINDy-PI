%% This file is used to test whether function "CalDerivative" could give us the correct result.
% Coded By: K
% Last Updated: 2019/06/19
%%
clc;clear;close all;

dt=0.001;T=5;
t=(0:dt:T)';

x=sin(t);
dx=cos(t);
ddx=-sin(t);
dddx=-cos(t);

dx_1=CalDerivative(x,dt,1);
ddx_1=CalDerivative(x,dt,2);
dddx_1=CalDerivative(x,dt,3);

%% Plot the result
figure(1)
hold on
plot(t,dx,'linewidth',3,'linestyle','-','color','black')
plot(t,dx_1,'linewidth',2.2,'linestyle','--','color','blue')
legend('True','Calculated')
grid on
title('dx')

figure(2)
hold on
plot(t,ddx,'linewidth',3,'linestyle','-','color','black')
plot(t,ddx_1,'linewidth',2.2,'linestyle','--','color','blue')
legend('True','Calculated')
grid on
title('ddx')


figure(3)
hold on
plot(t,dddx,'linewidth',3,'linestyle','-','color','black')
plot(t,dddx_1,'linewidth',2.2,'linestyle','--','color','blue')
legend('True','Calculated')
grid on
title('dddx')



