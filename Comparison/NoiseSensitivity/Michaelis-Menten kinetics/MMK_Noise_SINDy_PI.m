%% This file is the main file of camparing the noise sensitivity of iSINDy
% and SINDy-PI method. This file takes a long time to run.
%
% Last Update: 2019/07/11
% Coded By: K

%% Close all, clear all, clc
close all;clear all; clc;
addpath('./Functions')
addpath('./Datas')
set(0,'defaulttextInterpreter','latex')
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

% Define parameter for lamda
Sweep_Len=68;

% Determine whether you want to normalize your libary or not. 1 is yes and
% 0 is no.
NormalizeLib=0;

% Determine the model selection method
Model_Selection_Method=1;
% Determine the prediction step
Prediction_Steps=0;

% The follwoing parameters will determine the for loop and loop through the
% different noise level
d_percent=1;
percent=0;
percent_start=1;
% Must be smaller or equal to 24
percent_end=24;

% Determine a vector to store whether SINDy-PI identification is correct
is_Right=zeros(percent_end-percent_start+1,1);

% Turn off the rank deficient warning
MSGID='MATLAB:rankDeficientMatrix';
warning('off',MSGID);

% Determine how many iteration you need 
Final_pin=30;

% Choose which denoise method you want:
% 1 is adding Guassian noise to actual derivative, 2 is direct
% finite difference, 3 is TVRegDiff.
Noisy_Data_Generation_Method=3;

% Determine how many initial conditions you need
init_num=2400;

% Create the new directory to save the function files
FolderName=strcat('Result_DL_SINDy_Derivative_Method',num2str(Noisy_Data_Generation_Method),'_PridictionStep_',num2str(Prediction_Steps),'_IniNum_',num2str(init_num));
[fld_status, fld_msg, fld_msgID]=mkdir('TempFunctions_DL_SINDy');
[fld_status, fld_msg, fld_msgID]=mkdir(FolderName);
addpath('TempFunctions_DL_SINDy')

%% Now calculate the sensitivity of the iSINDy and SINDy-PI to the noise.

% Print the start process
fprintf('Start calculating ....\n\n\n')

for total=30:Final_pin
    percent=0;
    
    for percent_iter=percent_start:d_percent:percent_end
        % Set the pin one step forward
        percent=percent+1;
        
        % Set the noise level
        noise_level(percent,1)=Determine_Noise_Level(percent_iter);
        
        % Print the current noise level
        fprintf('Using the noise level as %i on SINDy-PI.\n',noise_level(percent,1)) 
        
        % Read files, load the training data and testing data
        File_Name=strcat('Data_Iter_',num2str(total),'_NoiseLevel_',num2str(noise_level(percent,1)),'_NoiseMethod_',num2str(Noisy_Data_Generation_Method),'_GoodLuck.mat');
        load(File_Name)
        [size1,size2]=size(xt_test);
        
        %% SINDy-PI
        fprintf('\n\n\n\n Using SINDy-PI now...\n')
        % First run the SINDy-PI
        for iter=1:n_state
            fprintf('\n \n Calculating the %i expression...\n',iter)
            
            % According to the previous parameter generate the left hand side guess
            [LHS_Data,LHS_Sym]=GuessLib(Data,dData(:,iter),u,Highest_Poly_Order_Guess,Highest_Trig_Order_Guess,Highest_U_Order_Guess);
            
            %Generate the corresponding data
            [SINDy_Data,SINDy_Struct]=SINDyLib(Data,dData(:,iter),u,Highest_Poly_Order,Highest_Trig_Order,Highest_U_Order,Highest_dPoly_Order);
            
            % Define the variable to store the ODE expression and score of each
            % expression. This is necessary to perform the parallel for loop.
            if percent==1 && iter==1
                ODEs_DL=cell(n_state,length(LHS_Sym),Sweep_Len);
                Score_DLS=NaN*ones(percent_end-percent_start+1,n_state,length(LHS_Sym),Sweep_Len);
                % Create a vector to store the result of whether current
                % ODE has the right structure
                Xi=cell(percent_end-percent_start+1,n_state,length(LHS_Sym),Sweep_Len);
            end
            
            % Run the for loop and try all the left hand guess
            for i=1:length(LHS_Sym)
                % Print the left hand side that we are testing
                fprintf('\t Testing the left hand side as %s:\n',char(LHS_Sym{i}))
                
                % Exclude the guess from SINDy library
                [RHS_Data,RHS_Struct]=ExcludeGuess(SINDy_Data,SINDy_Struct,LHS_Sym{i});
                
                % Sweep the values of lambda
                fprintf('\t\t Sweeping the values of lambda...\n')
                
                % Reassign the parameter for speed and memery
                LHS_Data_Dum=LHS_Data(:,i);
                LHS_Sym_Dum=LHS_Sym{i};
                dz_dum=dz(iter);
                
                % Swipe the parameters
                parfor j=1:Sweep_Len
                    % Generate the value of lambda
                    lambda=Get_Lambda(j);
                    
                    % Print the value of the lambda
                    fprintf('\t\t\t Testing lambda as %d...\n',lambda)
                    
                    % Perform the sparse regression problem
                    [Xi{percent,iter,i,j},ODE{j}]=sparsifyDynamics(RHS_Data,LHS_Data_Dum,LHS_Sym_Dum,lambda,N_iter,RHS_Struct,disp,NormalizeLib);
                    
                    % Store the result
                    try
                        % Try to store the ODE we find, perform sybolic calculation and solve for dX
                        digits(4)
                        ODEs_DL{iter,i,j}=vpa(solve(LHS_Sym_Dum==ODE{j},dz_dum));
                        
                        % Generate the ODE file
                        Generate_ODE_RHS(ODEs_DL{iter,i,j},n_state,n_control,strcat('TempFunctions_DL_SINDy/ParFor',num2str(j)));
                        % Calculate the accuracy of the file
                        Test_Score=0;
                        Test_Score=Get_Score(dxt_test,xt_test,u,Control,tspan,Shuffle,strcat('ParFor',num2str(j)),Prediction_Steps,dt,size1,size2,Model_Selection_Method);
                        Score_DLS(percent,iter,i,j)=Test_Score;
                    catch
                        ODEs_DL{iter,i,j}=NaN;
                        % If the ODEs_DL(percent,iter,i,j)=0, this means the solution does not
                        % exist. This represents the situation where dz=0.
                        Score_DLS(percent,iter,i,j)=NaN;
                    end
                end
            end
        end
    end
    %% Store the minumum of the SINDy-PI
    for pin1=1:size(Score_DLS,1)
        for pin3=1:size(Score_DLS,3)
            [min_DLS_Val_LHS(pin1,pin3),Index1(total,pin1,pin3)]=min(Score_DLS(pin1,1,pin3,:));
        end
        [min_DLS_Val(total,pin1),Index2(total,pin1)]=min(min_DLS_Val_LHS(pin1,:));
    end
    
    Coff_Error=zeros(Final_pin,size(Score_DLS,1));
    Dum3=zeros(9,1);
    Coffs1=[0.6;-3;0;0;0;-10/3;0;0;0];
    Coffs2=[0.18;-0.9;0;0;0;-0.3;0;0;0];
    for pin1=1:size(Score_DLS,1)
        Dum1=Index2(total,pin1);
        Dum2=Index1(total,pin1,Dum1);
        Dum3=Xi{pin1,1,Dum1,Dum2};
        Dum4=zeros(9,1);
        Dum4=Dum3~=0;
        if Dum1==1
            if Dum4==[1;1;0;0;0;1;0;0;0]
                is_Right(pin1,1)=is_Right(pin1,1)+1;
                Coff_Error(total,pin1)=0.3*norm(Coffs1-Dum3,1);
            else
                Coff_Error(total,pin1)=NaN;
            end
        elseif Dum1==2
            if Dum4==[1;1;0;0;0;1;0;0;0]
                is_Right(pin1,1)=is_Right(pin1,1)+1;
                Coff_Error(total,pin1)=norm(Coffs2-Dum3,1);
            else
                Coff_Error(total,pin1)=NaN;
            end
        elseif Dum1==3
            Coff_Error(total,pin1)=NaN;
        end
    end
    
    %% Plot the result for current iteration
    figure(1)
    plot(noise_level,min_DLS_Val(total,:),'o')
    drawnow
    grid on
    title('Peformance Comparison','FontSize',18)
    xlabel('Noise Level $\sigma$','FontSize', 18)
    ylabel('L2 Norm','FontSize', 18)
    set(gca,'FontSize',18);
    set(gcf,'Position',[100 100 600 400]);
    set(gcf,'PaperPositionMode','auto');
    set(gca, 'XScale', 'log')
    set(gca, 'YScale', 'log')
    %
    figure(2)
    clf
    plot(noise_level,mean(min_DLS_Val,1),'linewidth',2,'Color','green')
    drawnow
    grid on
    title('Peformance Comparison','FontSize',18)
    xlabel('Noise Level $\sigma$','FontSize', 18)
    ylabel('Average L2 Norm','FontSize', 18)
    set(gca,'FontSize',18);
    set(gcf,'Position',[100 100 600 400]);
    set(gcf,'PaperPositionMode','auto');
    set(gca, 'XScale', 'log')
    set(gca, 'YScale', 'log')
    %
    figure(3)
    plot(noise_level,is_Right,'o')
    drawnow
    grid on
    %legend('Success Rate')
    title('Success Rate','FontSize',18)
    xlabel('Noise Level $\sigma$','FontSize', 18)
    ylabel('Success Rate','FontSize', 18)
    ylim([0 total])
    set(gca,'FontSize',18);
    set(gcf,'Position',[100 100 600 400]);
    set(gcf,'PaperPositionMode','auto');
    set(gca, 'XScale', 'log')
    %%
    
    % Save the calculation result of current iteration
    cc=clock;
    ResultName=strcat(FolderName,'/DL_SINDY_Result_Iter_',num2str(total),'_Time_',num2str(total),'__',num2str(cc(3)),'_',num2str(cc(4)),'_',num2str(cc(5)),'_',num2str(round(cc(6))),'.mat');
    save(ResultName,'Score_DLS','Xi')
    
end




