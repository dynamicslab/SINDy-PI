%% This file is the main file of using the SINDy-PI method to
% infer the Single Pendulum on Cart Model.
%
% Date: 2019/07/26
% Coded By: K

%% Close all, clear all, clc
close all;clear all; clc;
dir
set(0,'defaulttextInterpreter','latex')
[status,message,messageid] = mkdir('Results');
addpath('Function');  
%% Simulate the budworm population growth and gather the simulation data

%Define the model parameters
M=1;m=1;L=1;g=9.81;

% Define noise level and add gaussian noise to the data
noise=0.02;

% Define whehter you have control, if you have it, please define it
Control=1;

%Define whether you want to shuffel the final data
Shuffle=0;

% Peform simulation
state0=[0.3;0;1;0];state0_test=[0.1;0;0.3;0];
dt=0.001;T=16;T_test=2;
tspan=0:dt:T;tspan_test=0:dt:T_test;
u=-0.2+0.5*sin(6*(tspan'));
u_test=-1+1*sin((tspan_test'))+3*sin(2*(tspan_test'));
[dData,Data]=Get_Sim_Data(@(t,y,u)SinglePendulum_ODE(t,y,u,M,m,L,g),state0,u,tspan,noise,Control,Shuffle);
dData(:,1)=Data(:,3);dData(:,2)=Data(:,4);
[dData_test,Data_test]=Get_Sim_Data(@(t,y,u)SinglePendulum_ODE(t,y,u,M,m,L,g),state0_test,u_test,tspan_test,noise,Control,Shuffle);
dData_test(:,1)=Data_test(:,3);dData_test(:,2)=Data_test(:,4);
%% Plot the simulation data
% close all
% HorizontalRange=5*L;
% VerticalRange=2.5*L;
% displayrate=1;
% cartoonOffLineAcce(Data(:,2),Data(:,1),L,HorizontalRange,VerticalRange,displayrate)
%%
% figure(1)
% plot(tspan,Data(:,1))
% title('Theta')
% grid on
% %
% figure(2)
% plot(tspan,Data(:,2))
% title('Position')
% grid on
% %
% figure(3)
% plot(tspan,Data(:,3))
% title('dTheta')
% grid on
% %
% figure(4)
% plot(tspan,Data(:,4))
% title('Velocity')
% grid on
% %
% figure(5)
% plot(tspan,u,'linewidth',3,'color','black')
% box('off')
% axis('on')

%% Now perform sparse regression of non-linear dynamics

% Get the number of states we have
[dtat_length,n_state]=size(Data);

% Define the control input(Should be zero in our example)
n_control=1;

% Choose whether you want to display actual ODE or not
disp_actual_ode=1;

% If the ODEs you want to display is the actual underlyting dynamics of the
% system, please set actual as 1
actual=1;

% Print the actual ODE we try to discover
digits(4)
z_vars=sym('z',[1,n_state]);
u_vars=sym('u',[1,n_control]);
d_vars=sym('dz',[1,n_state]);
ODEsP=SinglePendulum_ODE(0,z_vars,u_vars,M,m,L,g);
fprintf('The actual ODE of the system is/are :\n')
for i=1:n_state
    fprintf(strcat(char(d_vars(1,i)),'=',char(ODEsP(i,1)),'\n'));
end

% The implicit ODE has the following form:
% dz1=z3
% dz2=z4
% dz3=((981*sin(z1))/50 + u1*cos(z1) + z3^2*cos(z1)*sin(z1))/(cos(z1)^2 - 2)
% dz4=-(u1 + z3^2*sin(z1) + (981*cos(z1)*sin(z1))/100)/(cos(z1)^2 - 2)

% Create symbolic states
dz=sym('dz',[n_state,1]);
z=sym('z',[n_state,1]);
uc=sym('u',[n_control,1]);

% Now we first create the parameters of the function right hand side
Highest_Poly_Order_Guess=1;
Highest_Trig_Order_Guess=2;
Highest_U_Order_Guess=0;

% Then create the right hand side library parameters
Highest_Poly_Order=2;
Highest_Trig_Order=2;
Highest_U_Order=0;
Highest_dPoly_Order=1;

%% Define parameters for the sparese regression
lam=[1e-4;5e-4;1e-3;2e-3;3e-3;4e-3;5e-3;6e-3;7e-3;8e-3;9e-3;1e-2;2e-2;3e-2;4e-2;5e-2;...
    6e-2;7e-2;8e-2;9e-2;1e-1;2e-1;3e-1;4e-1;5e-1;6e-1;7e-1;8e-1;9e-1;1;1.5;2;2.5;3;3.5;4;4.5;5;...
    6;7;8;9;10;20;30;40;50;100;200];


N_iter=20;
disp=0;
LibPoly=[1;1;0;0];
LibTrig=[0;0;3;4];
%
GuessPoly=[1;1;0;0];
GuessTrig=[0;0;1;1];

% Normalize the library?
NormalizeLib=0;

% Define a cell matrix to store the variable
Xi_Final=cell(n_state,1);
Score_Final=cell(n_state,1);

for iter=1:n_state
    % Change the library base on what equation you want to work on
    Highest_Poly_Order=LibPoly(iter);
    Highest_Poly_Order_Guess=GuessPoly(iter);
    Highest_Trig_Order=LibTrig(iter);
    Highest_Trig_Order_Guess=GuessTrig(iter);
    
    fprintf('Calculating the %i expression \n',iter)
    
    % According to the previous parameter generate the left hand side guess
    [LHS_Data,LHS_Sym]=GuessLib(Data,dData(:,iter),iter,u,Highest_Poly_Order_Guess,Highest_Trig_Order_Guess,Highest_U_Order_Guess);
    
    %Generate the corresponding data
    [SINDy_Data,SINDy_Struct]=SINDyLib(Data,dData(:,iter),iter,u,Highest_Poly_Order,Highest_Trig_Order,Highest_U_Order,Highest_dPoly_Order);
    
    % Run the for loop and try all the left hand guess
    for i=1:length(LHS_Sym)
        fprintf('\n\t\t Testing the left hand side as %s... \n',cell2sym(LHS_Sym(1,i)))
        
        % Exclude the guess from SINDy-PI library
        [RHS_Data,RHS_Struct]=ExcludeGuess(SINDy_Data,SINDy_Struct,LHS_Sym{i});
        
        if iter==1
            Xi=cell(length(LHS_Sym),length(lam));
            Score=zeros(length(LHS_Sym),length(lam));
            ODE=cell(length(LHS_Sym),length(lam));
        end
        
        % Define dummy variables for parfor
        LHS_Sym_Dum=LHS_Sym{i};
        LHS_Data_Dum=LHS_Data(:,i);
        dz_Dum=dz(iter);
        dData_test_dum=dData_test(:,iter);
        for j=1:length(lam)
            fprintf('\n\t\t Testing the lambda as %d... \n',lam(j,1))
            % Perform the sparse regression problem
            [Xi{i,j},ODE{i,j}]=sparsifyDynamics(RHS_Data,LHS_Data_Dum,LHS_Sym_Dum,lam(j,1),N_iter,RHS_Struct,disp,NormalizeLib);
            
            % Perform sybolic calculation and solve for dX
            Eqn=LHS_Sym_Dum==ODE{i,j};
            digits(4)
            ODE_Guess=simplify(vpa(solve(Eqn,dz_Dum)));
            % Store the result of each valuation
            try
                ODEs(i,j)=ODE_Guess;
                func=matlabFunction(ODEs(i,j),'Vars',{z,uc});
                funcVal=func(Data_test',u_test')';
                Score(i,j)=norm(dData_test_dum-funcVal)/norm(dData_test(:,iter));
            catch
                ODEs(i,j)=0;
                Score(i,j)=NaN;
            end
        end
    end
    
    % Calculate the minimum score
    [minVal,minIndex]=min(Score,[],2);
    [minVal2,minIndex2]=min(minVal);
    
    % Store the best ODE approximation
    ODE_Best(iter,1)=ODEs(minIndex2,minIndex(minIndex2));
    
    % Print the result
    fprintf('\n \tThe SINDy-PI discovered Best ODE for %i equation is:\n',iter)
    digits(4)
    fprintf(strcat(char(dz(iter)),'=',char(simplify(ODE_Best(iter,1))),'\n'));
    
    % Save this result
    Xi_Final{iter,1}=Xi;
    Score_Final{iter,1}=Score;
    
    % Plot the score
    figure(iter)
    hold on
    for k=1:length(LHS_Sym)
        plot(lam,Score(k,:),'color',[1 iter*0.2 iter*0.25],'linewidth',2.5)
        scatter(lam,Score(k,:),100,'MarkerFaceColor',[1 iter*0.2 iter*0.25],'MarkerEdgeColor',[0 0 0],'linewidth',2.5)
        set(gca,'XScale','log')
    end
    
end


%% Now generate the ODE function file and test the accuracy of the
% identified system
Noise=0;

% Now generate this best guess ODE
Generate_ODE_RHS(ODE_Best(:,1),n_state,n_control);

% Simulate the system with the SINDy-PI identified ODE
%state0_test=[0.5;0.1;0.1;-0.1];
state0_test=[pi;0;0;0];
u_test=-0.5+0.2*sin((tspan_test'))+0.3*sin(2*(tspan_test'));
LHS_Sym=0;
[dData_test,Data_test]=Get_Sim_Data(@(t,y,u)SinglePendulum_ODE(t,y,u,M,m,L,g),state0_test,u_test,tspan_test,Noise,Control,Shuffle);
[d_Data_Es,Data_Es]=Get_Sim_Data(@(t,y,u)Sindy_ODE_RHS(t,y,u),state0_test,u_test,tspan_test,Noise,Control,Shuffle);

%% Print the Result
disp_best=1;
if disp_best==1
    fprintf('The SINDy-PI discovered Best ODE is:\n')
    digits(4)
    for i=1:n_state
        fprintf(strcat(char(dz(i)),'=',char(simplify(ODE_Best(i,1))),'\n'));
    end
end
%% Save the result
File_Name=strcat('Results/Noise_Level_',num2str(noise),'.mat');
save(File_Name,'Xi_Final','Score_Final','ODE_Best','state0_test',...
    'u_test','dData_test','Data_test','d_Data_Es','Data_Es','tspan_test')

%% Plot the simulation data
close all
% Create the new directory to save the plot
[fld_status, fld_msg, fld_msgID]=mkdir('Figures');
%
figure(1)
plot(tspan_test,Data_test(:,1),'linewidth',4.5,'Color','black')
hold on
plot(tspan_test,Data_Es(:,1),'linewidth',4.5,'linestyle','--','color','blue')
% legend('Actual Dynamics','Best Model')
% title('Validation','FontSize',24)
% xlabel('Time $(t)$','FontSize', 24)
% ylabel('$\theta(t)$','FontSize', 24)
set(gca,'FontSize',34);
grid on
h = gca; 
set(gca,'xticklabel',[])
set(gcf,'Position',[100 100 600 400]); 
set(gcf,'PaperPositionMode','auto');
box('on')
%print('-depsc2', '-loose', 'Figures/SinglePendulum_Thetat1.eps');

%
figure(2)
plot(tspan_test,Data_test(:,2),'linewidth',4.5,'Color','black')
hold on
plot(tspan_test,Data_Es(:,2),'linewidth',4.5,'linestyle','--','color','blue')
% legend('Actual Dynamics','Best Model')
% title('Validation','FontSize',18)
% xlabel('Time $(t)$','FontSize', 18)
% ylabel('$x(t)$','FontSize', 18)
set(gca,'FontSize',34);
grid on
%h = gca; h.XAxis.Visible = 'off';
set(gcf,'Position',[100 100 600 400]); 
set(gcf,'PaperPositionMode','auto');
box('on')
print('-depsc2', '-loose', 'Figures/SinglePendulum_x.eps');

%%
figure(3)
plot(tspan_test,Data_test(:,3),'linewidth',2,'Color','green')
hold on
plot(tspan_test,Data_Es(:,3),'linewidth',2,'linestyle','--','color','blue')
legend('Actual Dynamics','Best Model')
title('Validation','FontSize',18)
xlabel('Time $(t)$','FontSize', 18)
ylabel('$\dot{\theta}(t)$','FontSize', 18)
set(gca,'FontSize',18);
grid on

set(gcf,'Position',[100 100 600 400]); 
set(gcf,'PaperPositionMode','auto');
print('-depsc2', '-loose', 'Figures/SinglePendulum_dThetat.eps');

%%
figure(4)
plot(tspan_test,Data_test(:,4),'linewidth',2,'Color','green')
hold on
plot(tspan_test,Data_Es(:,4),'linewidth',2,'linestyle','--','color','blue')
legend('Actual Dynamics','Best Model')
title('Validation','FontSize',18)
xlabel('Time $(t)$','FontSize', 18)
ylabel('$\dot{x}(t)$','FontSize', 18)
set(gca,'FontSize',18);
grid on

set(gcf,'Position',[100 100 600 400]); 
set(gcf,'PaperPositionMode','auto');
print('-depsc2', '-loose', 'Figures/SinglePendulum_dx.eps');

%
figure(5)
plot(Data_test(:,1),Data_test(:,3),'linewidth',2,'Color','green')
hold on
plot(Data_Es(:,1),Data_Es(:,3),'linewidth',2,'linestyle','--','color','blue')
legend('Actual Dynamics','Best Model')
title('Validation','FontSize',18)
xlabel('$\theta(t)$','FontSize', 18)
ylabel('$\dot{\theta}(t)$','FontSize', 18)
set(gca,'FontSize',18);
grid on

set(gcf,'Position',[100 100 600 400]); 
set(gcf,'PaperPositionMode','auto');
print('-depsc2', '-loose', 'Figures/SinglePendulum_theta_vs_dtheta.eps');

%
figure(6)
plot(Data_test(:,2),Data_test(:,4),'linewidth',2,'Color','green')
hold on
plot(Data_Es(:,2),Data_Es(:,4),'linewidth',2,'linestyle','--','color','blue')
legend('Actual Dynamics','Best Model')
title('Validation','FontSize',18)
xlabel('$x(t)$','FontSize', 18)
ylabel('$\dot{x}(t)$','FontSize', 18)
set(gca,'FontSize',18);
grid on

set(gcf,'Position',[100 100 600 400]); 
set(gcf,'PaperPositionMode','auto');
print('-depsc2', '-loose', 'Figures/SinglePendulum_x_vs_dx.eps');










