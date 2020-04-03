%% This file is the main file of using the SINDy-PI method to
% infer the ODE of the double pendulum.
%
% Date: 2019/05/06
% Coded By: K

%% Close all, clear all, clc
close all;clear all; clc;
[status,msg] = mkdir('Results');
addpath('Functions')
set(0,'defaulttextInterpreter','latex')
%% Simulate the double pendulum and gather the simulation data

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
k1=0;k2=0;

%Define the simulation time length
Tf=10;dt=0.001;tspan=0:dt:Tf;
T_test=3;tspan_test=0:dt:T_test;

%Define the inital state of the pendulum: theta1, theta2, dtheta1, dtheta2
state0=[pi+1.2;pi-0.6;0;0];
state0_test=[pi-1;pi-0.4;0.3;0.4];

% Define noise level and add gaussian noise to the data
noise=0;

% Define whehter you have control, if you have it, please define it
Control=0;u=0;

%Define whether you want to shuffel the final data
Shuffle=0;

% Run the ODE files and gather the simulation data
[dData,Data]=Get_Sim_Data(@(t,y)DouPenODE(t,y,m1,m2,a1,a2,L,I1,I2,k1,k2,g),state0,u,tspan,noise,Control,Shuffle);

[dData_test,Data_test]=Get_Sim_Data(@(t,y)DouPenODE(t,y,m1,m2,a1,a2,L,I1,I2,k1,k2,g),state0_test,u,tspan,noise,Control,Shuffle);

%% Plot the data
figure(1)
plot(tspan,Data(:,3),'linewidth',3,'color','black')
box('off')
axis('on')
set(gca,'FontSize',24)

figure(2)
plot(tspan,Data(:,4),'linewidth',3,'color','black')
box('off')
axis('on')
set(gca,'FontSize',24)

% figure(3)
% plot(Data(:,3),dData(:,3),'linewidth',3,'color','black')
% box('off')
% axis('off')
% 
% figure(4)
% plot(Data(:,4),dData(:,4),'linewidth',3,'color','black')
% box('off')
% axis('off')

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
Print_ODEs(@(t,y)DouPenODE(t,y,m1,m2,a1,a2,L,I1,I2,k1,k2,g),n_state,n_control,disp_actual_ode,actual);

% Create symbolic states
dz=sym('dz',[n_state,1]);

% Now we first create the parameters of the function right hand side
Highest_Poly_Order_Guess=1;
Highest_Trig_Order_Guess=1;
Highest_U_Order_Guess=0;

% Then create the right hand side library parameters
Highest_Poly_Order=1;
Highest_Trig_Order=4;
Highest_U_Order=0;
Highest_dPoly_Order=1;

%% Define parameters for the sparese regression
lam=[1e-4;5e-4;1e-3;2e-3;3e-3;4e-3;5e-3;6e-3;7e-3;8e-3;9e-3;1e-2;2e-2;3e-2;4e-2;5e-2;...
    6e-2;7e-2;8e-2;9e-2;1e-1;2e-1;3e-1;4e-1;5e-1;6e-1;7e-1;8e-1;9e-1;1;1.5;2;2.5;3;3.5;4;4.5;5;...
    6;7;8;9;10;20;30;40;50;100;200];

N_iter=20;
disp=0;
NormalizeLib=0;

for iter=1:n_state
    fprintf('\n \n Calculating the %i expression...\n',iter)
    
    % According to the previous parameter generate the left hand side guess
    [LHS_Data,LHS_Sym]=GuessLib(Data,dData(:,iter),iter,u,Highest_Poly_Order_Guess,Highest_Trig_Order_Guess,Highest_U_Order_Guess);
    
    %Generate the corresponding data
    [SINDy_Data,SINDy_Struct]=SINDyLib(Data,dData(:,iter),iter,u,Highest_Poly_Order,Highest_Trig_Order,Highest_U_Order,Highest_dPoly_Order);
    
    % Run the for loop and try all the left hand guess
    for i=1:length(LHS_Sym)
        if iter==1 && i==1
            Xi=cell(n_state,length(LHS_Sym),length(lam));
            ODE=cell(n_state,length(LHS_Sym),length(lam));
            ODEs=cell(n_state,length(LHS_Sym),length(lam));
        end

        % Print the left hand side that we are testing
        fprintf('\t Testing the left hand side as %s:\n',char(LHS_Sym{i}))
        
        % Exclude the guess from SINDy library
        [RHS_Data,RHS_Struct]=ExcludeGuess(SINDy_Data,SINDy_Struct,LHS_Sym{i});
        
        parfor j=1:length(lam)
            % Select the sparse threashold
            lambda=lam(j);
    
            % Perform the sparse regression problem
            [Xi{iter,i,j},ODE{iter,i,j}]=sparsifyDynamics(RHS_Data,LHS_Data(:,i),LHS_Sym{i},lambda,N_iter,RHS_Struct,disp,NormalizeLib);
            
            % Perform sybolic calculation and solve for dX
            digits(6)
            ODE_Guess=vpa(solve(LHS_Sym{i}==ODE{iter,i,j},dz(iter)));
            
            % Print the discovered ODE
            fprintf(strcat('\t The corresponding ODE we found is: ',char(dz(iter,1)),'=',char((ODE_Guess)),'\n \n'));
            
            % Store the result
            ODEs{iter,i,j}=ODE_Guess;
        end
    end
end

%% Now generate the ODE function file and test the accuracy of the
% identified system

fprintf('\v Start calculating the best model that could represent the training data...\n \n')
for iter=1:n_state
    % Print which expression are you working on
    fprintf('\t Calculating the best model for the %d expression...\n',iter)
    
    for i=1:length(LHS_Sym)
        % Print the process
        fprintf('\t Calculating the score of previously found ODE on the test data, %d %% finished. \n',round((i/length(LHS_Sym))*100))
        
        for j=1:length(lam)
            % If the previous ODE is 0, set the score as NaN, else calculate
            % it.
            if isempty(ODEs{iter,i,j})
                ODE_Not_Exist=1;
                Score(iter,i,j)=NaN;
            else
                % Generate the ODE file
                Generate_ODE_RHS(ODEs{iter,i,j},n_state,n_control);
                % Calculate the accuracy of the file
                Score(iter,i,j)=Get_Score(dData_test(:,iter),Data_test,u,Control,tspan_test,state0_test,Shuffle);
            end
        end
        
        % Get the best lambda
        [minVal1(iter,i),minIndex1(iter,i)]=min(Score(iter,i,:));
        
    end
    
    % Get the best score and use this ODE file
    [minVal2(iter,1),minIndex2(iter,1)]=min(minVal1(iter,:));
    
    % Store the best ODE
    ODE_Best(iter,1)=ODEs{iter,minIndex2(iter,1),minIndex1(iter,minIndex2(iter,1))};
    
    % Print the Result
    fprintf('\n\n\n\t The SINDy-PI discovered Best ODE for the %d expression is:\n',iter)
    fprintf('\t %s = %s \n\n\n',char(dz(iter)),char((ODE_Best(iter,1)))')
end

%% Now generate this best guess ODE and print its result
disp_best=1;
if disp_best==1
    fprintf('\n\n\n\v The SINDy-PI discovered Best ODE for the whole system is:\n')
    digits(4)
    for iter=1:n_state
        fprintf(strcat('\v ******\v\t',char(dz(iter)),'=',char(simplify(ODE_Best(iter,1))),'\n'));
    end
end

% Also Print the actual ODE for comparison
digits(4)
fprintf('\n\n\n')
Print_ODEs(@(t,y)DouPenODE(t,y,m1,m2,a1,a2,L,I1,I2,k1,k2,g),n_state,n_control,disp_actual_ode,actual);

%% Get the simulation result
% Generate the ODE file
fprintf('\n\n\n\v Generating the Best Model for comparision...\n')
Generate_ODE_RHS(ODE_Best(:,1),n_state,n_control);

% Define a new test data for the comparison
Noise_test=0;
state0_test=[pi+0.3;pi-0.5;0;0];
[dData_Es,Data_Es]=Get_Sim_Data(@(t,z)Sindy_ODE_RHS(t,z,u),state0_test,u,tspan_test,Noise_test,Control,Shuffle);
[dData_test,Data_test]=Get_Sim_Data(@(t,y)DouPenODE(t,y,m1,m2,a1,a2,L,I1,I2,k1,k2,g),state0_test,u,tspan_test,Noise_test,Control,Shuffle);

%% Save the result
File_Name=strcat('Results/DoublePendulum_NoiseLevel_',num2str(noise),'.mat');
save(File_Name,'Score','ODEs','ODE_Best','state0_test','Data_Es','dData_Es',...
    'Data_test','dData_test','Xi','tspan_test')

%% Plot the simulation data
% Show process
fprintf('\n Simulation finished, plotting the result...... \n')

close all
% Create the new directory to save the plot
[fld_status, fld_msg, fld_msgID]=mkdir('Figures');

%%
close all
figure(1)
plot(tspan_test,Data_test(:,1),'linewidth',4.5,'Color','black')
hold on
plot(tspan_test,Data_Es(:,1),'linewidth',4.5,'linestyle','--','color','blue')
% title('Validation')
% xlabel('Time $(t)$')
% ylabel('$\theta_1$')
% legend('Actual Dynamics','Approximated Dynamics')
set(gca,'XTickLabel',[])
set(gca,'FontSize',34);
grid on

set(gcf,'Position',[100 100 600 400]);
set(gcf,'PaperPositionMode','auto');
%print('-depsc2', '-loose', 'Figures/DoublePendulum_Theta1.eps');

%
figure(2)
plot(tspan_test,Data_test(:,2),'linewidth',4.5,'Color','black')
hold on
plot(tspan_test,Data_Es(:,2),'linewidth',4.5,'linestyle','--','color','blue')
% title('Validation')
% xlabel('Time $(t)$')
% ylabel('$\theta_2$')
% legend('Actual Dynamics','Approximated Dynamics')
set(gca,'FontSize',34);
grid on

set(gcf,'Position',[100 100 600 400]);
set(gcf,'PaperPositionMode','auto');
print('-depsc2', '-loose', 'Figures/DoublePendulum_Theta2.eps');

%%
figure(3)
plot(Data_test(:,1),dData_test(:,1),'linewidth',4.5,'Color','green')
hold on
plot(Data_Es(:,1),dData_Es(:,1),'linewidth',4.5,'linestyle','--','color','blue')
% title('Phase Plot: $\theta_1$ vs $\dot{\theta_1}$')
% xlabel('$\theta_1$')
% ylabel('$\dot{\theta_1}$')
% legend('Actual Dynamics','Approximated Dynamics')
set(gca,'FontSize',24);
grid on

set(gcf,'Position',[100 100 600 400]);
set(gcf,'PaperPositionMode','auto');
print('-depsc2', '-loose', 'Figures/DoublePendulum_PhasePlot_Theta1_vs_Theta2.eps');

%
figure(4)
plot(Data_test(:,2),dData_test(:,2),'linewidth',4.5,'Color','green')
hold on
plot(Data_Es(:,2),dData_Es(:,2),'linewidth',4.5,'linestyle','--','color','blue')
% title('Phase Plot: $\theta_2$ vs $\dot{\theta_2}$')
% xlabel('$\theta_2$')
% ylabel('$\dot{\theta_2}$')
% legend('Actual Dynamics','Approximated Dynamics')
set(gca,'FontSize',24);
grid on

set(gcf,'Position',[100 100 600 400]);
set(gcf,'PaperPositionMode','auto');
print('-depsc2', '-loose', 'Figures/DoublePendulum_PhasePlot_dTheta1_vs_dTheta2.eps');

