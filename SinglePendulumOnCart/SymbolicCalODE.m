%% This file will symbolic solve for the ODE of the Single Pendulum on a Cart system.
% Coded By: K
% Data: 2019/04/27
%%
clc;clear;close all
% Define the symbolic symbols you need
syms x0(t) dx0(t) ddx0 theta(t) dtheta(t) ddtheta M m L F k0 k1 g I1 dum
% Set some variables to zero to simplify the model
I1=0; k0=0; k1=0;
% First define the position of the pendulum mass center
x1=x0+L*sin(theta)
y1=L*cos(theta)

% Now define the velocity of the cart and moving mass
v0=diff(x0,t)
v1x=diff(x1,t)
v1y=diff(y1,t)

% Substitute the variables
v0=subs(v0,diff(x0,t),dx0);
v1x=subs(v1x,diff(x0,t),dx0);
v1x=subs(v1x,diff(theta,t),dtheta);
v1y=subs(v1y,diff(x0,t),dx0);
v1y=subs(v1y,diff(theta,t),dtheta);

simplify(v1x^2+v1y^2)
simplify(v1x)
simplify(v1y)

% Now define the kinetic energy of the system
EngK=simplify(0.5*M*v0^2+0.5*m*simplify(v1x^2+v1y^2)+0.5*I1*diff(theta,t)^2)
EngP=m*g*(cos(theta))*L

% Define the Rayley disspation function
Damp=0.5*k0*v0^2+0.5*k1*dtheta^2;

% Now define the Lagrangian
Lag=simplify(EngK-EngP)

%% Now get the second order ODE
the1Part1=diff(subs(Lag,dtheta,dum),dum);
the1Part1=subs(the1Part1,dum,dtheta);
the1Part1=diff(the1Part1,t);
%
the1Part2=diff(subs(Lag,theta,dum),dum);
the1Part2=subs(the1Part2,dum,theta);
% 
the1Part3=diff(subs(Damp,dtheta,dum),dum);
the1Part3=subs(the1Part3,dum,dtheta);
eqn1=the1Part1-the1Part2+the1Part3
%
eqn1=subs(eqn1,diff(dx0(t), t),ddx0);
eqn1=subs(eqn1,diff(x0(t), t, t),ddx0);
eqn1=subs(eqn1,diff(theta(t), t),dtheta);
eqn1=subs(eqn1,diff(dtheta(t), t),ddtheta);

% Second equation
the2Part1=diff(subs(Lag,dx0,dum),dum);
the2Part1=subs(the2Part1,dum,dx0);
the2Part1=diff(the2Part1,t);
%
the2Part2=diff(subs(Lag,x0,dum),dum);
the2Part2=subs(the2Part2,dum,x0);
% 
the2Part3=diff(subs(Damp,dx0,dum),dum);
the2Part3=subs(the2Part3,dum,dx0);
eqn2=the2Part1-the2Part2+the2Part3-F
%
eqn2=subs(eqn2,diff(dx0(t), t),ddx0);
eqn2=subs(eqn2,diff(x0(t), t, t),ddx0);
eqn2=subs(eqn2,diff(theta(t), t),dtheta);
eqn2=subs(eqn2,diff(dtheta(t), t),ddtheta);

% Final eqn1 and eqn2
Eqn1=simplify(eqn1)==0 
Eqn2=simplify(eqn2)==0

%% Solve for ddtheta and ddx0
eqns=[Eqn1 Eqn2];
vars=[ddtheta ddx0];
[Ans1 Ans2]=solve(eqns,vars);

Ans1=simplify(Ans1)
Ans2=simplify(Ans2)

%% Now rewrite the second order system into first order system
syms z1 z2 z3 z4
thisisddtheta1=simplify(subs(Ans1,[theta(t) x0(t) dtheta(t) dx0(t)],[z1 z2 z3 z4]))
thisisddtheta2=simplify(subs(Ans2,[theta(t) x0(t) dtheta(t) dx0(t)],[z1 z2 z3 z4]))

%% Save this file into Matlab Function
% dtheta=z3, dx=z4, ddtheta=thisistheta1, ddx=thisistheta2
Eqns=[z3;z4;thisisddtheta1;thisisddtheta2]

matlabFunction(Eqns,'File','SinglePendulum_ODE','Vars',{t,[z1 z2 z3 z4],F,M,m,L,g},'Optimize',true)







