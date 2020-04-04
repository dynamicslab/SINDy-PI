%% This file is the main file of using the SINDy-PI method to
% infer the Yeast Glycolysis Model. We will figure out what is the minimum
% data length needed for SINDy-PI to accurately discover the six-th state
% of the Yeast Glycolysis Model. This file will be used to swipe through
% different data length.
%
% Date: 2019/06/13
% Coded By: K

%% Close all, clear all, clc
close all;clear all; clc;
set(0,'defaulttextInterpreter','latex')
addpath('Functions')
addpath('Datas')
%% Define some parameters
% Define whehter you have control, if you have it, please define it
n_control=0;u=0;

% Run the ODE files and gather the simulation data.
% (We use the same data for both the iSINDy method and SINDy-PI method for better comparision)
% (Baseline data for the data length comparison of state 1,2,3,4,5,6,7)
% The data is noise clean. It is generated using 900 different initial
% conditions with 51 points of each initial condition.
load('TrainingData.mat')

% Choose whether you want to display actual ODE or not
disp_actual_ode=1;

% If the ODEs you want to display is the actual underlyting dynamics of the
% system, please set actual as 1
actual=1;

% Define how many states we have in our example
n_state=7;

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

% Set the parameter normLib=1 to normalize the librry
normLib=1;

% Set how amny iterations you want for the sparse regression
N_iter=10;

% Set whether you want to display the ODE or not
disp=0;

% Set the library size
Highest_Poly_Order_Lib=[6;6;3;3;3;6;3];

% Determine how many percent of data you want.
percent_start=1;d_percent=0.005;percent_end=12;


% Determine which states you want to discover
Which_State=6;

% Determine the left hand side guess number
if Which_State==1 || Which_State==2 || Which_State==6
    LHS_Num=2;
else
    LHS_Num=1;
end

% Get the poly-order
Highest_Poly_Order=Highest_Poly_Order_Lib(Which_State);

% Determine which LHS you want to use, the first one or the second one
LHS_Pin=1;
%%
tic
% Create the library data using all the data points
[SINDy_Data_Full,SINDy_Struct]=SINDyLib(xt,dxt(:,Which_State),Which_State,u,Highest_Poly_Order,Highest_Trig_Order,Highest_U_Order,Highest_dPoly_Order);
fprintf('\n\t Original library creation finished, using %i seconds...\n',toc)

% Create left hand side guess
[LHS_Data_Full,LHS_Sym]=GuessLib(xt,dxt(:,Which_State),Which_State,u,Highest_Poly_Order_Guess,Highest_Trig_Order_Guess,Highest_U_Order_Guess);

tic
% Create the right hand side, exclude the guess from SINDy library
[RHS_Data_Full,RHS_Struct]=ExcludeGuess(SINDy_Data_Full,SINDy_Struct,LHS_Sym{LHS_Pin});
LHS_Data_Full_Dum=LHS_Data_Full(:,LHS_Pin);
fprintf('\t Library without LHS guess creation finished, using %i seconds...\n',toc)

% Determine sparsity parameter
lambda=0.1;

% Create the new directory to save the result
FolderName=strcat('Result_DL_SINDy_Data_Length_Compare_State_',num2str(Which_State),'_LHS_Guess_',num2str(LHS_Pin),'_SwipeLength_',num2str(lambda));
[fld_status, fld_msg, fld_msgID]=mkdir(FolderName);

% Set up dummy variables for parfor
LHS_Sym_Dum=LHS_Sym{LHS_Pin};
dz_dum=dz(Which_State);
%% determine the percentage you want to use
if LHS_Pin==1
    Percent_Iter=[0.30 0.31 0.32 0.33 0.34 0.345 0.35 0.355 0.36 0.365 0.37 0.38 0.39 0.40];
else
    Percent_Iter=[0.02 0.03 0.04 0.05 0.055 0.0575 0.06 0.0625 0.065 0.0675 0.07 0.0725 0.075 0.08 0.09 0.1 0.2];
end

for per=1:size(Percent_Iter,2)
    percent=Percent_Iter(per);

    % Start!
    for Total_Run=1:10
        fprintf('\n \n Get the result for the %i time...\n',Total_Run)
        
        % Print the state we are testing
        fprintf('\n \t Calculating the %i expression...\n',Which_State)
        
        % Print the left hand side that we are testing
        fprintf('\t Testing the left hand side as %s:\n',char(LHS_Sym_Dum))
        
        % Define the new data length
        new_length=round((percent)*length(xt));
        
        % Shuffel the original data
        Sequence=randperm(size(SINDy_Data_Full,1));
        Sequence_Trimed=Sequence(1:round((percent)*length(xt)));
        
        RHS_Data=RHS_Data_Full(Sequence_Trimed,:);
        LHS_Data=LHS_Data_Full_Dum(Sequence_Trimed,1);
        
        fprintf('\n \t Data preparation finished, using %i seconds...\n',toc)
        
        % We fix the lambda and change the data usage.
        tic
        fprintf('\n \t Testing the percentage as %i ...\n',(percent*100))
        
        % Perform the sparse regression problem
        Xi=sparsifyDynamics_simplified(RHS_Data,LHS_Data,lambda,N_iter,normLib);
        fprintf('\t Uses %i seconds...\n',toc)
        
        % Save the calculation result of current iteration
        fprintf('\n\t Saving the result...\n')
        tic
        cc=clock;
        ResultName=strcat(FolderName,'/DL_SINDY_Data_Length_',num2str(Total_Run),'__','LHS',num2str(LHS_Pin),'_',num2str(cc(3)),'_',num2str(cc(4)),'_',num2str(cc(5)),'_',num2str(round(cc(6))),'_P3','.mat');
        save(ResultName,'Xi','lambda','percent')
        fprintf('\n\t Saving finished! Using %i seconds...\n',toc)
    end
end





