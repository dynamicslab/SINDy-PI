% This file will use SINDy-PI on MMK equation.
%
% Date: 2019/07/19
% Coded By: K

%% Close all, clear all, clc
close all;clear all; clc;
addpath('./Functions')
set(0,'defaulttextInterpreter','latex')
%% Simulate the budworm population growth and gather the simulation data

% Define the system parameters
jx=0.6;Vmax=1.5;Km=0.3;

% Determine the simulation time step and time span
dt=0.1; T=5; tspan=0:dt:T;

% Define whehter you have control, if you have it, please define it
Control=0;u=0;

%Define whether you want to shuffel the final data
Shuffle=0;

% Load the clean data
load('MMK_Simulation.mat')
[size1,size2]=size(xt_test);
Test_Data=reshape(xt_test,[],1);
Test_dData=reshape(dxt_test,[],1);

%% Set up some parameters
% Get the number of states we have
n_state=1;

% Define the control input(Should be zero in our example)
n_control=0;

% Choose whether you want to display actual ODE or not
disp_actual_ode=1;

% If the ODEs you want to display is the actual underlyting dynamics of the
% system, please set actual as 1
actual=1;

% Print the actual ODE we try to discover
Print_ODEs(@(t,y)MMK_ODE(t,y,jx,Vmax,Km),n_state,n_control,disp_actual_ode,actual);

% Create symbolic states
dz=sym('dz',[n_state,1]);

% Now we first create the parameters of the function right hand side
Highest_Poly_Order_Guess=1;
Highest_Trig_Order_Guess=1;
Highest_U_Order_Guess=0;

% Then create the right hand side library parameters
Highest_Poly_Order=4;
Highest_Trig_Order=0;
Highest_U_Order=0;
Highest_dPoly_Order=1;

% Define parameters for the sparese regression
N_iter=20;
disp=0;

%% Now calculate the sensitivity of the iSINDy and SINDy-PI to the noise.

% Create the new directory to save the function files
[fld_status, fld_msg, fld_msgID]=mkdir('TempFunctions');
[fld_status, fld_msg, fld_msgID]=mkdir('Result_DL_SINDy_Small_Noise');
addpath('TempFunctions')

% Print the start process
fprintf('Start calculating ....\n\n\n')

% Define parameters
Sweep_Len=30;
NormalizeLib=0;

% Set the noise level
noise_level=0;

% Print the current noise level
fprintf('Using the noise level as %i on SINDy-PI.\n',noise_level)

% Assign the value to the new variables
init_num=2400;
Noisy_Data_Generation_Method=3;
Gap=200;
[Data,dData,iter_length]=Generate_Noisy_Data(noise_level,init_num,dt,T,Noisy_Data_Generation_Method,Gap);
Prediction_Steps=0;
Model_Selection_Method=1;

%% SINDy-PI
%tic
fprintf('\n\n\n\n Using SINDy-PI now...\n')
% First run the SINDy-PI

for iter=1:n_state
    fprintf('\n \n Calculating the %i expression...\n',iter)
    
    % According to the previous parameter generate the left hand side guess
    [LHS_Data,LHS_Sym]=GuessLib(Data,dData(:,iter),u,Highest_Poly_Order_Guess,Highest_Trig_Order_Guess,Highest_U_Order_Guess);
    
    %Generate the corresponding data
    [SINDy_Data,SINDy_Struct]=SINDyLib(Data,dData(:,iter),u,Highest_Poly_Order,Highest_Trig_Order,Highest_U_Order,Highest_dPoly_Order);
    
    % Run the for loop and try all the left hand guess
    for i=1:length(LHS_Sym)
        % Print the left hand side that we are testing
        fprintf('\t Testing the left hand side as %s:\n',char(LHS_Sym{i}))
        
        % Exclude the guess from SINDy library
        [RHS_Data,RHS_Struct]=ExcludeGuess(SINDy_Data,SINDy_Struct,LHS_Sym{i});
        
        % Sweep the values of lambda
        fprintf('\t\t Sweeping the values of lambda...\n')
        
        if iter==1 && i==1
            Xi=cell(length(LHS_Sym),Sweep_Len);
            ODE=cell(length(LHS_Sym),Sweep_Len);
            ODEs_DL=cell(length(LHS_Sym),Sweep_Len);
            Score_DLS=zeros(n_state,length(LHS_Sym),Sweep_Len);
        end
        
        LHS_Data_Dum=LHS_Data(:,i);
        LHS_Sym_Dum=LHS_Sym{i};
        dz_dum=dz(iter);
        parfor j=1:Sweep_Len
            % Generate the value of lambda
            lambda=Get_Lambda(j);
            
            % Print the value of the lambda
            fprintf('\t\t\t Testing lambda as %d...\n',lambda)
            
            % Perform the sparse regression problem
            [Xi{i,j},ODE{i,j}]=sparsifyDynamics(RHS_Data,LHS_Data_Dum,LHS_Sym_Dum,lambda,N_iter,RHS_Struct,disp,NormalizeLib);
            
            % Store the result
            try
                % Try to store the ODE we find, perform sybolic calculation and solve for dX
                digits(4)
                ODEs_DL{i,j}=vpa(solve(LHS_Sym_Dum==ODE{i,j},dz_dum));
                
                % Generate the ODE file
                Generate_ODE_RHS(ODEs_DL{i,j},n_state,n_control,strcat('TempFunctions/ParFor',num2str(j)));
                % Calculate the accuracy of the file
                Test_Score=0;
                Test_Score=Get_Score(Test_dData,Test_Data,u,Control,tspan,Shuffle,strcat('ParFor',num2str(j)),Prediction_Steps,dt,size1,size2,Model_Selection_Method);
                Score_DLS(iter,i,j)=Test_Score;
            catch
                ODEs_DL{i,j}=NaN;
                % If the ODEs_DL(percent,iter,i,j)=0, this means the solution does not
                % exist. This represents the situation where dz=0.
                Score_DLS(iter,i,j)=NaN;
            end
        end
    end
end
%%
Right=0;
for i=1:Sweep_Len
    Dum3=Xi{1,i};
    Dum4=Dum3~=0;
    if sum(abs(Dum4-[1;1;0;0;0;1;0;0;0]))==0
        Right=Right+1;
        Dum3;
        ODEs_DL{1,i};
    end
    Dum3=Xi{2,i};
    Dum4=Dum3~=0;
    if sum(abs(Dum4-[1;1;0;0;0;1;0;0;0]))==0
        Right=Right+1;
        Dum3;
        ODEs_DL{2,i};
    end
end

%  Get the best model
[minVal1,Index1]=min(Score_DLS(1,1,:));
[minVal2,Index2]=min(Score_DLS(1,2,:));
[minVal3,Index3]=min(Score_DLS(1,3,:));

Best_Index=[Index1 Index2 Index3];

[minVal,Index]=min([minVal1 minVal2 minVal3]);

%% Print the ODE
fprintf('\n\n The best model discovered by SINDy-PI is:\n %s\n',ODEs_DL{Index,Best_Index(Index)})

%% Plot
for j=1:Sweep_Len
    x_axis(j,1)=Get_Lambda(j);
end

close all
figure(1)
hold on
% ss1=scatter(x_axis,squeeze(Score_DLS(1,1,:)),100,'filled',...
%     'MarkerEdgeColor',[0 0 0],...
%     'MarkerFaceColor',1*[1 1 1],...
%     'LineWidth',1.5)
% ss1.Marker = 'o';

ss2=scatter(x_axis,squeeze(Score_DLS(1,2,:)),100,'filled',...
    'MarkerEdgeColor',[0 0 0],...
    'MarkerFaceColor',1*[0 0.4 1],...
    'LineWidth',1.5)
ss2.Marker = 'o';

ss3=scatter(x_axis,squeeze(Score_DLS(1,3,:)),100,'filled',...
    'MarkerEdgeColor',[0 0 0],...
    'MarkerFaceColor',1*[1 0.4 0],...
    'LineWidth',1.5)
ss3.Marker = 'o';

box('on')
set(gca,'FontSize',24)
set(gca,'XScale','log','YScale','log')
grid on



