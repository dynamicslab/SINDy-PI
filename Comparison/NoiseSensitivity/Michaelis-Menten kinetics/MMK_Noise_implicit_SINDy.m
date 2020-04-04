%% This file is the main file of camparing the noise sensitivity of iSINDy
% and SINDy-PI method. This file takes a long time to run.
%
% Last Update: 2019/07/17
% Coded By: K

%% Close all, clear all, clc
close all;clear all; clc;
addpath('./Functions')
addpath('Datas')
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
Highest_Trig_Order_Guess=0;
Highest_U_Order_Guess=0;

% Then create the right hand side library parameters
Highest_Poly_Order=4;
Highest_Trig_Order=0;
Highest_U_Order=0;
Highest_dPoly_Order=1;

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

% Define a matrix to store the null space vector
Xi_ns=cell(percent_end-percent_start+1,1);

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
FolderName=strcat('Result_i_SINDy_Derivative_Method',num2str(Noisy_Data_Generation_Method),'_PridictionStep_',num2str(Prediction_Steps),'_IniNum_',num2str(init_num));
[fld_status, fld_msg, fld_msgID]=mkdir('TempFunctions_iSINDy');
[fld_status, fld_msg, fld_msgID]=mkdir(FolderName);
addpath('TempFunctions_iSINDy')


%% Now calculate the sensitivity of the iSINDy and SINDy-PI to the noise.

% Print the start process
fprintf('Start calculating ....\n\n\n')

for total=1:Final_pin
    percent=0;
    
    for percent_iter=percent_start:d_percent:percent_end
        % Set the pin one step forward
        percent=percent+1;
        
        % Set the noise level
        noise_level(percent,1)=Determine_Noise_Level(percent_iter);
        
        % Print the current noise level
        fprintf('Using the noise level as %i on implicit-SINDy.\n',noise_level(percent,1))
        
        % Read files, load the training data and testing data
        File_Name=strcat('Data_Iter_',num2str(total),'_NoiseLevel_',num2str(noise_level(percent,1)),'_NoiseMethod_',num2str(Noisy_Data_Generation_Method),'_GoodLuck.mat');
        load(File_Name)
        [size1,size2]=size(xt_test);
        
        %% Run implicit-SINDy (some function in this section is obtained from the github code iSINDy )
        fprintf('\n\n\n\n USing implicit-SINDy now...\n')
        
        % Sweep through the states
        %tic
        for iter=1:n_state
            fprintf('\n \n Calculating the %i expression...\n',iter)
            
            % Build library data
            [iSINDy_Data,iSINDy_Struct]=SINDyLib(Data,dData(:,iter),u,Highest_Poly_Order,Highest_Trig_Order,Highest_U_Order,Highest_dPoly_Order);
            
            % Use Null-Space method and set the tolerance
            tol = 1e-5;
            
            % Plotting option
            pflag=0;
            
            % Print process
            fprintf('\n \n \t Calculating the null space...\n')
            try
                % Sparse Regression
                [Xi_ns{percent,1}, indTheta, lambdavec, numterms, errorv] = ADMpareto(iSINDy_Data, tol, pflag);
                Xi_dum=zeros(10,28);
                Xi_dum=Xi_ns{percent,1};
                %Define the variable to store the ODE expression and score of each
                % expression. This is necessary to perform the parallel for loop.
                ODEs_iS=cell(n_state,size(Xi_dum,2)-1);
                if percent==1 && iter==1
                    Score_iS=NaN*ones(percent_end-percent_start+1,size(Xi_dum,2)-1);
                end
                dz_dum=dz(iter);
                
                % Sweep the null space and store the score
                parfor j=1:(size(Xi_dum,2)-1)
                    
                    % Print which vector we are working on
                    fprintf('\t\t Sweeping the %i null sapce vector...\n',j)
                    
                    % Validate each sparse vector
                    try
                        % Solve for the ODE expression
                        ODEs_iS{iter,j}=solve(vpa(cell2sym(iSINDy_Struct)*Xi_dum(:,j))==0,dz_dum);
                        
                        % Store the previous result into a matrix
                        Generate_ODE_RHS(ODEs_iS{iter,j},n_state,n_control,strcat('TempFunctions_iSINDy/ParForiS',num2str(j)));
                        
                        % Calculate the accuracy of the file
                        Test_Score=0;
                        Test_Score=Get_Score(dxt_test,xt_test,u,Control,tspan,Shuffle,strcat('ParForiS',num2str(j)),Prediction_Steps,dt,size1,size2,Model_Selection_Method);
                        Score_iS(percent,j)=Test_Score;
                    catch
                        ODEs_iS{iter,j}=NaN;
                        Score_iS(percent,j)=NaN;
                    end
                end
            catch
                fprintf('\n\t\t\tSomething is wrong....\n')
            end
            
        end
        
    end
    
    % Store the result for this iteration
    
    % Store the minimum of implicit SINDy
    for pin1=1:size(Score_iS,1)
        [min_iS_Val(total,pin1),Index1(total,pin1)]=min(Score_iS(pin1,:));
        
    end
    
    % Test whether the correct structure is identified
    for pin1=1:size(Score_iS,1)
        Xi_dum=Xi_ns{pin1,1};
        Coffs=Xi_dum(:,Index1(total,pin1));
        Coffs_dum=Coffs~=0;
        if Coffs_dum==[1;1;0;0;0;1;1;0;0;0]
            is_Right(pin1,1)=is_Right(pin1,1)+1;
            Coffs_norm=Coffs/Coffs(7);
            Coffs_error(total,pin1)=norm(Coffs_norm-[-0.18;0.9;0;0;0;0.3;1;0;0;0]);
        elseif Coffs_dum==[0;1;1;0;0;0;1;1;0;0]
            is_Right(pin1,1)=is_Right(pin1,1)+1;
            Coffs_norm=Coffs/Coffs(8);
            Coffs_error(total,pin1)=norm(Coffs_norm-[0;-0.18;0.9;0;0;0;0.3;1;0;0]);
        elseif Coffs_dum==[0;0;1;1;0;0;0;1;1;0]
            is_Right(pin1,1)=is_Right(pin1,1)+1;
            Coffs_norm=Coffs/Coffs(9);
            Coffs_error(total,pin1)=norm(Coffs_norm-[0;0;-0.18;0.9;0;0;0;0.3;1;0]);
        elseif Coffs_dum==[0;0;0;1;1;0;0;0;1;1]
            is_Right(pin1,1)=is_Right(pin1,1)+1;
            Coffs_norm=Coffs/Coffs(10);
            Coffs_error(total,pin1)=norm(Coffs_norm-[0;0;0;-0.18;0.9;0;0;0;0.3;1]);
        end
    end
    
    for pin1=1:size(Score_iS,1)
        try
            Num=sum(Coffs_error(:,pin1)~=0);
            Ave_Coff_Error(pin1,1)=sum(Coffs_error(:,pin1))/Num;
        catch
            Ave_Coff_Error(pin1,1)=NaN;
        end
    end
    
    %% Plot the result
    figure(1)
    plot(noise_level,log(1+mean(min_iS_Val,1)),'o','linewidth',2,'color','blue')
    drawnow
    grid on
    title('Peformance Comparison','FontSize',18)
    xlabel('Noise Level $\sigma$','FontSize', 18)
    ylabel('Average L2 Norm','FontSize', 18)
    set(gca,'FontSize',18);
    set(gcf,'Position',[100 100 600 400]);
    set(gcf,'PaperPositionMode','auto');
    
    figure(2)
    plot(noise_level,is_Right/total,'o')
    drawnow
    grid on
    title('Success Rate','FontSize',18)
    xlabel('Noise Level $\sigma$','FontSize', 18)
    ylabel('Success Rate','FontSize', 18)
    ylim([0 1])
    set(gca,'FontSize',18);
    set(gcf,'Position',[100 100 600 400]);
    set(gcf,'PaperPositionMode','auto');
    set(gca, 'XScale', 'log')
    
    figure(3)
    plot(noise_level,Ave_Coff_Error,'o')
    drawnow
    grid on
    title('Average Parameter Error','FontSize',18)
    xlabel('Noise Level $\sigma$','FontSize', 18)
    ylabel('L2 norm','FontSize', 18)
    ylim([0 total])
    set(gca,'FontSize',18);
    set(gcf,'Position',[100 100 600 400]);
    set(gcf,'PaperPositionMode','auto');
    set(gca, 'XScale', 'log')
    
    % Save the calculation result of current iteration
    cc=clock;
    ResultName=strcat(FolderName,'/implicit_SINDY_Result',num2str(total),'__',num2str(cc(3)),'_',num2str(cc(4)),'_',num2str(cc(5)),'_',num2str(round(cc(6))),'.mat');
    save(ResultName,'Score_iS','Xi_ns')
    
end




