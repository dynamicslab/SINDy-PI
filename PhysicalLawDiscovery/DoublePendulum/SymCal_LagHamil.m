%% This file will calculate the Euler-Lagrange's equation of the double pendulum 
% and it's Hamiltonian.
%
% Coded By: K
% Last Updated: 05/06/2019
%% Clear all, close all
clc;clear;

%% Define the parameters and variables that we need
syms M m1 m2 L1 L2 I1 I2 k0 k1 k2 g x(t) the1(t) the2(t) a1 a2 dum
syms dx(t) dthe1(t) dthe2(t) ddx(t) ddthe1(t) ddthe2(t) dumddthe1 dumddthe2

%k1=0;k2=0;

%% Calculate the intermediate variables 
% Calculate the position of mass center of each pendulum arm
x0=x;
x1=x0+a1*sin(the1);
x2=x0+L1*sin(the1)+a2*sin(the2);

y1=a1*cos(the1);
y2=L1*cos(the1)+a2*cos(the2);

% Calculate the velocity of the mass center
v0=diff(x0,t);
v1x=diff(x1,t);
v1y=diff(y1,t);
v2x=diff(x2,t);
v2y=diff(y2,t);

%Substitute variables
v0=subs(v0,diff(x(t), t),dx);
v1x=subs(v1x,diff(x(t), t),dx);
v1y=subs(v1y,diff(x(t), t),dx);
v2x=subs(v2x,diff(x(t), t),dx);
v2y=subs(v2y,diff(x(t), t),dx);

%Subs dthe1
v1x=subs(v1x,diff(the1(t), t),dthe1);
v1y=subs(v1y,diff(the1(t), t),dthe1);
v2x=subs(v2x,diff(the1(t), t),dthe1);
v2y=subs(v2y,diff(the1(t), t),dthe1);

%Subs dthe2
v1x=subs(v1x,diff(the2(t), t),dthe2);
v1y=subs(v1y,diff(the2(t), t),dthe2);
v2x=subs(v2x,diff(the2(t), t),dthe2);
v2y=subs(v2y,diff(the2(t), t),dthe2);

%% Calculate  the Lagrangian and Hamiltonian
%Define the Lagrangian and Damping coeificient
Damp=0.5*k0*v0^2+0.5*k1*dthe1^2+0.5*k2*(dthe2-dthe1)^2;
Lag=0.5*M*v0^2+0.5*m1*(v1x^2+v1y^2)+0.5*m2*(v2x^2+v2y^2)+0.5*(I1)*dthe1^2+0.5*I2*dthe2^2-m1*g*y1-m2*g*y2;
LagEng=0.5*M*v0^2+0.5*m1*(v1x^2+v1y^2)+0.5*m2*(v2x^2+v2y^2)+0.5*(I1)*dthe1^2+0.5*I2*dthe2^2+m1*g*y1+m2*g*y2;

%% Substitute the intermediate variables
% Subsitute the first one
the1Part1=diff(subs(Lag,dthe1,dum),dum);
the1Part1=subs(the1Part1,dum,dthe1);
the1Part1=diff(the1Part1,t);
%
the1Part2=diff(subs(Lag,the1,dum),dum);
the1Part2=subs(the1Part2,dum,the1);
% 
the1Part3=diff(subs(Damp,dthe1,dum),dum);
the1Part3=subs(the1Part3,dum,dthe1);
eqn1=the1Part1-the1Part2+the1Part3;

% Subsitute the second one
the2Part1=diff(subs(Lag,dthe2,dum),dum);
the2Part1=subs(the2Part1,dum,dthe2);
the2Part1=diff(the2Part1,t);
%
the2Part2=diff(subs(Lag,the2,dum),dum);
the2Part2=subs(the2Part2,dum,the2);
% 
the2Part3=diff(subs(Damp,dthe2,dum),dum);
the2Part3=subs(the2Part3,dum,dthe2);
eqn2=the2Part1-the2Part2+the2Part3;

% Now substitute the second derivative to variables
eqn1=subs(eqn1,diff(dx(t), t),ddx);
eqn1=subs(eqn1,diff(x(t), t, t),ddx);

eqn1=subs(eqn1,diff(the1(t), t),dthe1);
eqn1=subs(eqn1,diff(the2(t), t),dthe2);

eqn1=subs(eqn1,diff(dthe1(t), t),ddthe1);
eqn1=subs(eqn1,diff(dthe2(t), t),ddthe2);

eqn2=subs(eqn2,diff(dx(t), t),ddx);
eqn2=subs(eqn2,diff(x(t), t, t),ddx);

eqn2=subs(eqn2,diff(the1(t), t),dthe1);
eqn2=subs(eqn2,diff(the2(t), t),dthe2);

eqn2=subs(eqn2,diff(dthe1(t), t),ddthe1);
eqn2=subs(eqn2,diff(dthe2(t), t),ddthe2);

%Nowsolve for ddthe1 and ddthe2
eqn1=subs(eqn1,ddthe1,dumddthe1);
eqn1=subs(eqn1,ddthe2,dumddthe2);
%
eqn2=subs(eqn2,ddthe1,dumddthe1);
eqn2=subs(eqn2,ddthe2,dumddthe2);
%
eqn1=(simplify(eqn1)==0);
eqn2=(simplify(eqn2)==0);
%
eqns=[eqn1 eqn2];
Lag1=simplify(eqn1);
Lag2=simplify(eqn2);

%% Prrtify the expression
% Define the new variable: theta1=z1, theta2=z2, dtheta1=dz1, dtheta2=dz2,
% ddtheta1 =ddz1, ddtheta2=ddz2. For the mounted double pendulum we do not 
% consider the cart movement and set variables related to the x position to 
% zeros.  
syms z1 z2 dz1 dz2 ddz1 ddz2

% Prretify and calculate the total energy (Hamiltonian).
Lag1=subs(Lag1,[the1(t) the2(t) dthe1(t) dthe2(t) dumddthe1 dumddthe2 ddx(t)],[z1 z2 dz1 dz2 ddz1 ddz2 0]);
Lag2=subs(Lag2,[the1(t) the2(t) dthe1(t) dthe2(t) dumddthe1 dumddthe2 ddx(t)],[z1 z2 dz1 dz2 ddz1 ddz2 0]);
Hamil=subs(LagEng,[the1(t) the2(t) dthe1(t) dthe2(t) dumddthe1 dumddthe2 dx(t) ddx(t)],[z1 z2 dz1 dz2 ddz1 ddz2 0 0]);

% Print the EOM
fprintf('\n\n\n **********The symbolic form of the expression are:**********\n\n\n')
fprintf('\v Print the EOM of the double pendulum \n')
fprintf('\n\t ***The first EOM is:***\n')
fprintf('\n\t %s \n',char(simplify(Lag1)))

fprintf('\n\t ***The second EOM is:***\n')
fprintf('\n\t %s \n',char(simplify(Lag2)))

% Print the Hamiltonian
fprintf('\n\t ***The Hamiltonian of the system is:***\n')
fprintf('\n\t %s \n',char(simplify(Hamil)))

%% We calculate the system ODE in this sub section
% The ODE of the system is calcultaed by solving the EOM.
vars=[dumddthe1 dumddthe2];
[Ans1 Ans2]=solve(eqns,vars);

%ddthe1
Ans1=simplify(Ans1);

%ddthe2
Ans2=simplify(Ans2);

% Pretify the result
thisisddtheta1=simplify(subs(Ans1,[the1(t) the2(t) dthe1(t) dthe2(t) ddx(t)],[z1 z2 dz1 dz2 0]));
thisisddtheta2=simplify(subs(Ans2,[the1(t) the2(t) dthe1(t) dthe2(t) ddx(t)],[z1 z2 dz1 dz2 0]));

% Print the result of ODE
fprintf('\v The ODE of the double pendulum \n')
fprintf('\n\t ***The first ODE is:***\n')
fprintf('\n\t %s \n',char(simplify(thisisddtheta1)))

fprintf('\n\t ***The second ODE is:***\n')
fprintf('\n\t %s \n',char(simplify(thisisddtheta2)))

%% Take in the actual value of the parameters
%Load the system parameters, those parameters are based on the actual system
load("EstimatedValueDou.mat")

% Substitude the variable as number
Lag1Num=subs(Lag1,[m1 m2 a1 a2 I1 I2 k1 k2 g L1],[parsEsDou(1) parsEsDou(2) parsEsDou(3) parsEsDou(4) parsEsDou(5) parsEsDou(6) parsEsDou(7) parsEsDou(8) 9.81 0.2667]);
Lag2Num=subs(Lag2,[m1 m2 a1 a2 I1 I2 k1 k2 g L1],[parsEsDou(1) parsEsDou(2) parsEsDou(3) parsEsDou(4) parsEsDou(5) parsEsDou(6) parsEsDou(7) parsEsDou(8) 9.81 0.2667]);
HamilNum=subs(Hamil,[m1 m2 a1 a2 I1 I2 k1 k2 g L1],[parsEsDou(1) parsEsDou(2) parsEsDou(3) parsEsDou(4) parsEsDou(5) parsEsDou(6) parsEsDou(7) parsEsDou(8) 9.81 0.2667]);

thisisddtheta1Num=simplify(subs(thisisddtheta1,[m1 m2 a1 a2 I1 I2 k1 k2 g L1],[parsEsDou(1) parsEsDou(2) parsEsDou(3) parsEsDou(4) parsEsDou(5) parsEsDou(6) parsEsDou(7) parsEsDou(8) 9.81 0.2667]));
thisisddtheta2Num=simplify(subs(thisisddtheta2,[m1 m2 a1 a2 I1 I2 k1 k2 g L1],[parsEsDou(1) parsEsDou(2) parsEsDou(3) parsEsDou(4) parsEsDou(5) parsEsDou(6) parsEsDou(7) parsEsDou(8) 9.81 0.2667]));

%% Use this when the k1 and k2 is zero

% Lag1Num=subs(Lag1,[m1 m2 a1 a2 I1 I2 k1 k2 g L1],[parsEsDou(1) parsEsDou(2) parsEsDou(3) parsEsDou(4) parsEsDou(5) parsEsDou(6) 0 0 9.81 0.2667]);
% Lag2Num=subs(Lag2,[m1 m2 a1 a2 I1 I2 k1 k2 g L1],[parsEsDou(1) parsEsDou(2) parsEsDou(3) parsEsDou(4) parsEsDou(5) parsEsDou(6) 0 0 9.81 0.2667]);
% HamilNum=subs(Hamil,[m1 m2 a1 a2 I1 I2 k1 k2 g L1],[parsEsDou(1) parsEsDou(2) parsEsDou(3) parsEsDou(4) parsEsDou(5) parsEsDou(6) 0 0 9.81 0.2667]);
% 
% thisisddtheta1Num=simplify(subs(thisisddtheta1,[m1 m2 a1 a2 I1 I2 k1 k2 g L1],[parsEsDou(1) parsEsDou(2) parsEsDou(3) parsEsDou(4) parsEsDou(5) parsEsDou(6) 0 0 9.81 0.2667]));
% thisisddtheta2Num=simplify(subs(thisisddtheta2,[m1 m2 a1 a2 I1 I2 k1 k2 g L1],[parsEsDou(1) parsEsDou(2) parsEsDou(3) parsEsDou(4) parsEsDou(5) parsEsDou(6) 0 0 9.81 0.2667]));

%%
% Print the EOM
digits(4)
fprintf('\n\n\n **********The numerical form of the expression are:**********\n\n\n')
fprintf('\n\v Print the EOM of the double pendulum \n')
fprintf('\n\t ***The first EOM is:***\n')
fprintf('\n\t %s \n',char(simplify(vpa(Lag1Num/0.02751))))
%
fprintf('\n\t The second EOM is:***\n')
fprintf('\n\t %s \n',char(simplify(vpa(Lag2Num/0.006544))))

% Print the Hamiltonian
fprintf('\n\v The Hamiltonian of the double pendulum \n')
fprintf('\n\t ***The Hamiltonian of the system is:***\n')
fprintf('\n\t %s \n',char(simplify(vpa(HamilNum))))

% Print the result of ODE
fprintf('\n\v The ODE of the double pendulum \n')
fprintf('\n\t ***The first ODE is:***\n')
fprintf('\n\t %s \n',char(simplify(vpa(thisisddtheta1Num))))

fprintf('\n\t ***The second ODE is:***\n')
fprintf('\n\t %s \n',char(simplify(vpa(thisisddtheta2Num))))


