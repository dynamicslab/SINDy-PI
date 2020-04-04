%% This file will generate the necessary simulation data needed for the selection of
% SINDy model. In short, it will generate the test data. Note that we will
% run this file only once. The test data will remains the same during the
% model selection process.
%
% Date: 2019/06/03
% Coded By: K

%% Close all, clear all, clc
close all;clear all; clc;

%% Simulate the budworm population growth and gather the simulation data

% Define the system parameters
jx=0.6;Vmax=1.5;Km=0.3;

% Determine the simulation time step and time span
dt=0.1; T=5; tspan=0:dt:T;

% Define noise level and add gaussian noise to the data
noise=0;

% Define whehter you have control, if you have it, please define it
Control=0;u=0;

%Define whether you want to shuffel the final data
Shuffle=0;


%% Generate the test data
% Define how many initial condisions you need
init_num_test=100;

% Now run a for loop and store the simulation data
dxt_test=[];xt_test=[];
for iter=1:init_num_test
    x0_test(iter,1)=abs(((round(iter/10))+1)*rand);
    [dxt_dum,xt_dum]=Get_Sim_Data(@(t,y)MMK_ODE(t,y,jx,Vmax,Km),x0_test(iter,1),u,tspan,noise,Control,Shuffle);
    dxt_test=[dxt_test,dxt_dum];
    xt_test=[xt_test,xt_dum];
end
%% Save the result
save('TestData.mat','dxt_test','xt_test','x0_test','dt','T','tspan')


