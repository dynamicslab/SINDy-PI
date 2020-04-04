%% This file is the main file of using the SINDy-PI method to
% infer the Yeast Glycolysis Model.
%
% Date: 2019/05/06
% Coded By: K

%% Close all, clear all, clc
close all;clear all; clc;
set(0,'defaulttextInterpreter','latex')
addpath('Functions')
addpath('Datas')
%% Simulate the budworm population growth and gather the simulation data
%Define the simulation time length
Tf=5;dt=0.1;tspan=0:dt:Tf;
Tspan=0:dt:Tf;
T_test=5;
tspan_test=0:dt:T_test;

% Define noise level and add gaussian noise to the data
noise=0;
Noise_test=0;

% Define whehter you have control, if you have it, please define it
Control=0;u=0;

%Define whether you want to shuffel the final data
Shuffle=0;

% Run the ODE files and gather the simulation data(we use the same data generated for the iSINDy method for better comparision)
% (Baseline data for the data length comparison of state 1,2,3,4,5,6,7)
load('TrainingData.mat')

% This one is prepared for the 6th state, it might need this much data for the iSINDy algorithm. 
% (Baseline data for the data length comparison of state 6)
%load('TrainingData_Large.mat')

% Define the new data length
percent=0.1;new_length=round(percent*length(xt));

% Shuffel the original data
Sequence=randperm(size(xt,1));
xt=xt(Sequence,:);
dxt=dxt(Sequence,:);

% Assign the value to the new variables
Data=xt(1:new_length,:);dData=dxt(1:new_length,:);

%% Now perform sparse regression of non-linear dynamics

% Get the number of states we have
[dtat_length,n_state]=size(Data);

% Define the control input(Should be zero in our example)
n_control=0;

% Choose whether you want to display actual ODE or not
disp_actual_ode=1;

% If the ODEs you want to display is the actual underlyting dynamics of the
% system, please set actual as 1
actual=1;

% Print the actual ODE we try to discover
digits(4)
Print_ODEs(@(t,y)YeastGlycolysis_ODE(t,y),n_state,n_control,disp_actual_ode,actual);

% Create symbolic states
dz=sym('dz',[n_state,1]);

% Now we first create the parameters of the function right hand side
Highest_Poly_Order_Guess=0;
Highest_Trig_Order_Guess=0;
Highest_U_Order_Guess=0;

% Then create the right hand side library parameters
Highest_Trig_Order=0;
Highest_U_Order=0;
Highest_dPoly_Order=1;

%% Define parameters for the sparese regression

% Set the parameter normLib=1 to normalize the librry
normLib=1;

% You could manually set the threashold here and test the result of
% SINDy-PI.
lam=[0.5;0.5;0.6;0.8;0.4;0.2;0.2];
N_iter=20;
disp=0;
Highest_Poly_Order_Lib=[6;6;3;3;3;6;3];

for iter=3:3
    Highest_Poly_Order=Highest_Poly_Order_Lib(iter);
    
    fprintf('\n \n Calculating the %i expression...\n',iter)
    
    % Set the threashold
    lambda=lam(iter);
    
    % According to the previous parameter generate the left hand side guess
    [LHS_Data,LHS_Sym]=GuessLib(Data,dData(:,iter),iter,u,Highest_Poly_Order_Guess,Highest_Trig_Order_Guess,Highest_U_Order_Guess);
    
    %Generate the corresponding data
    [SINDy_Data,SINDy_Struct]=SINDyLib(Data,dData(:,iter),iter,u,Highest_Poly_Order,Highest_Trig_Order,Highest_U_Order,Highest_dPoly_Order);
    %%
    % Run the for loop and try all the left hand guess
    for i=1:length(LHS_Sym)
        % Print the left hand side that we are testing
        fprintf('\t Testing the left hand side as %s:\n',char(LHS_Sym{i}))
        
        % Exclude the guess from SINDy library
        [RHS_Data,RHS_Struct]=ExcludeGuess(SINDy_Data,SINDy_Struct,LHS_Sym{i});
        
        % Perform the sparse regression problem
        [Xi,ODE]=sparsifyDynamics(RHS_Data,LHS_Data(:,i),LHS_Sym{i},lambda,N_iter,RHS_Struct,disp,normLib);
        
        % Perform sybolic calculation and solve for dX
        Eqn=LHS_Sym{i}==ODE;
        digits(6)
        ODE_Guess=vpa(solve(Eqn,dz(iter)));
        
        % Print the discovered ODE
        fprintf(strcat('\t The corresponding ODE we found is: ',char(dz(iter,1)),'=',char(simplify(ODE_Guess)),'\n \n'));
        
        % Store the result
        ODEs(iter,i)=ODE_Guess;
    end
end

