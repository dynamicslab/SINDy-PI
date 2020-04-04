%% This file is the ODE file used to simulate the Michaelis-Menten kinetics.
% Date: 04/24/2019
% Coded By: K

% Here Vmax is rhe maximum rate of reaction time
% Km is the concentration of half-maximal reaction rate
% 

%% ODE functions
function dx=MMK_ODE(t,x,jx,Vmax,Km)
dx=jx-(Vmax*x(1,:)./(Km+x(1,:)));






