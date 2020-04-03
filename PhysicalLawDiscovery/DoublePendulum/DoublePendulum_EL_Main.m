%% This file is the main file of using the SINDy-PI method to
% infer the physical law of given data. We will try to model the
% Hamiltonian of the double pendulum in this file.
%
% Date: 2019/05/07
% Coded By: K

%% Close all, clear all, clc
close all;clear all; clc;
set(0,'defaulttextInterpreter','latex')
addpath('Functions')
%% Simulate the budworm population growth and gather the simulation data

%Load the system parameters, those parameters are based on the actual system
load("EstimatedValueDou.mat")

%Pendulum mass
m1=parsEsDou(1);m2=parsEsDou(2);

%Pendulum center of mass
a1=parsEsDou(3);a2=parsEsDou(4);

%Pendulum arm inertial
I1=parsEsDou(5);I2=parsEsDou(6);

%Gravity constant and the first pendulum arm length
g=9.81;L=0.2667;

%Dapming ratio
k1=parsEsDou(7);k2=parsEsDou(8);
%k1=0;k2=0;

%Define the simulation time length
Tf=15;dt=0.001;tspan=0:dt:Tf;
T_test=3;tspan_test=0:dt:T_test;

%Define the inital state of the pendulum: theta1, theta2, dtheta1, dtheta2
state0=[pi-0.6;pi-0.4;0;0];
state0_test=[pi-0.2;pi-0.1;0.1;0.1];

% Define noise level and add gaussian noise to the data
%noise=0.05;
noise=0;

% Define whehter you have control, if you have it, please define it
Control=0;u=0;

%Define whether you want to shuffel the final data
Shuffle=0;

% Run the ODE files and gather the simulation data
[dData,Data]=Get_Sim_Data(@(t,y)DouPenODE(t,y,m1,m2,a1,a2,L,I1,I2,k1,k2,g),state0,u,tspan,noise,Control,Shuffle);

[dData_test,Data_test]=Get_Sim_Data(@(t,y)DouPenODE(t,y,m1,m2,a1,a2,L,I1,I2,k1,k2,g),state0_test,u,tspan,noise,Control,Shuffle);

% Get the number of states we have
[dtat_length,n_state]=size(Data);

% Define the control input(Should be zero in our example)
n_control=0;

% Create symbolic states
dz=sym('dz',[n_state,1]);
z=sym('z',[n_state,1]);
syms z1 z2 z3 z4 dz1 dz2 dz3 dz4

% Define expression of Hamiltonian 
HamilExp=1.045*cos(z1) + 0.3268*cos(z2) + 0.01376*dz1^2 + 0.003272*dz2^2 + 0.008885*dz1*dz2*cos(z1 - z2);

% Save the function file of the Hamiltonian
matlabFunction(HamilExp,'File','DoublePendulumEnergy','Optimize',true,'Vars',{[z1,z2,z3,z4],[dz1,dz2,dz3,dz4]});

% Calculate the system energy of training data
SysEng=DoublePendulumEnergy(Data,dData);

%% Plot the data
close all
figure(1)
plot(tspan,Data(:,1),'linewidth',2.5,'color','black')
set(gca,'FontSize',24)
box('off')
%grid on

figure(2)
plot(tspan,Data(:,2),'linewidth',2.5,'color','black')
set(gca,'FontSize',24)
box('off')
%grid on

figure(3)
plot(tspan,Data(:,3),'linewidth',2.5,'color','black')
set(gca,'FontSize',24)
box('off')
%grid on

figure(4)
plot(tspan,Data(:,4),'linewidth',2.5,'color','black')
set(gca,'FontSize',24)
box('off')
%grid on
% 
% figure(5)
% plot(tspan,SysEng)
% grid on

%% Now perform sparse regression of non-linear dynamics
Highest_Poly_Order_Guess=0;
Highest_Trig_Order_Guess=0;
Highest_U_Order_Guess=0;

% Then create the right hand side library parameters
Highest_Poly_Order=0;
Highest_Trig_Order=0;
Highest_U_Order=0;
Highest_dPoly_Order=0;

%% Define parameters for the sparese regression
N_iter=20;
disp=0;
NormalizeLib=0;

for iter=1:1
    fprintf('\n \n Calculating the %i expression...\n',iter)
    
    % Select the sparse threashold
    lambda=0.002;
    
    % According to the previous parameter generate the left hand side guess
    [LHS_Data,LHS_Sym]=GuessLib(Data,dData,iter,u,Highest_Poly_Order_Guess,Highest_Trig_Order_Guess,Highest_U_Order_Guess);
    
    %Generate the corresponding data
    [SINDy_Data,SINDy_Struct]=SINDyLib(Data,dData,iter,u,Highest_Poly_Order,Highest_Trig_Order,Highest_U_Order,Highest_dPoly_Order);
    
    % Run the for loop and try all the left hand guess
    for i=1:length(LHS_Sym)
        % Print the left hand side that we are testing
        fprintf('\t Testing the left hand side as %s:\n',char(LHS_Sym{i}))
        
        % Exclude the guess from SINDy library
        [RHS_Data,RHS_Struct]=ExcludeGuess(SINDy_Data,SINDy_Struct,LHS_Sym{i});
        
        % Perform the sparse regression problem
        [Xi,ODE]=sparsifyDynamics(RHS_Data,LHS_Data(:,i),LHS_Sym{i},lambda,N_iter,RHS_Struct,disp,NormalizeLib);
        
        % Print the discovered ODE
        digits(4)
        fprintf(strcat('\t The corresponding equation we found is: ',char(vpa(ODE)),'  =   ',char(LHS_Sym{i}),'\n \n'));
        
        % Store the result
        ODEs(iter,i)=ODE;
        
    end
end

