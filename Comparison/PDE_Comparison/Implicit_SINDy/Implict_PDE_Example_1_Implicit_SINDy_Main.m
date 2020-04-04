%% This file is used to demonstrate that implicit-SINDy could be a bad choice to discover the implicit PDE.
% Coded By: K
% Last Updated: 2019/07/15
%% Clear all
clc;clear;close all

%% Create folder and read function files
addpath('Functions')
[fld_status, fld_msg, fld_msgID]=mkdir('TempFunctions');
addpath('TempFunctions')

%% We first generate the PDE data
% Define the time span
dt=0.01;T=20;
tspan=0:dt:T;

% Define the spatial domain
L=50;

% Define the discretization point
n=256;

% Whether you want to plot it or not
PDE_Plot=1;

Compare=0;

% Run simulation and gather data
[utsol,usol,x,k]=Example_1(dt,T,L,n,PDE_Plot,Compare);
dx=x(2)-x(1);

%% Now prepare the library data

% Prepare ux,uxx
ux=[];uxx=[];uxxx=[];
ux_dum=[];uxx_dum=[];uxxx_dum=[];

% Take derivative using finite difference
for kk=1:size(usol,2)
    ux(kk,:)=CalDerivative(real(usol(:,kk)),dx,1);
    uxx(kk,:)=CalDerivative(real(usol(:,kk)),dx,2);
    uxxx(kk,:)=CalDerivative(real(usol(:,kk)),dx,3);
end

% Reshape into a vector
Ux=reshape(ux',[],1);
Uxx=reshape(uxx',[],1);
Uxxx=reshape(uxxx',[],1);

%% Prepare the ut
ut=[];utt=[];
for xpos=1:size(usol,1)
    dummy=real(usol(xpos,:));
    ut(:,xpos)=CalDerivative(dummy,dt,1);
    utt(:,xpos)=CalDerivative(dummy,dt,2);
end

% Reshape
Ut=reshape(ut',[],1);
Utt=reshape(utt',[],1);

%% Prepare the u
u=real(usol)';

% Reshape
U=reshape(u',[],1);

%% Plot the result
close all
if PDE_Plot==1
    figure(1)
    %
    surf(x,tspan(1:1600),real(usol(:,1:1600)'),'FaceColor',[0.5 0.6 0.9],'EdgeColor','none','FaceLighting','none')
    %shading interp
    %colormap(gray)
    hold on
    %
    surf(x,tspan(1600:end),real(usol(:,1600:end)'),'FaceColor',[0.9 0.5 0.6],'EdgeColor','none','FaceLighting','none')
    %shading interp
    %
    Spacer1=5;
    Spacer2=50;
    %
    ss=surf(x(1:Spacer1:size(x,2)),tspan(1:Spacer2:size(tspan,2)),real(usol(1:Spacer1:size(x,2),1:Spacer2:size(tspan,2))'),'FaceAlpha',0)
    ss.EdgeColor='black'
    %
    set(gca,'FontSize',18);
    set(gcf,'Position',[100 100 600 400]);
    set(gcf,'PaperPositionMode','auto');
    grid on
    
    view(42,35)
end

%% Here, we prepare our data into testing set and training set
Data_Length=size(U,1);
Random_Sequence=randperm(Data_Length);

% Permute the original data
U=U(Random_Sequence);Ut=Ut(Random_Sequence);Utt=Utt(Random_Sequence);
Ux=Ux(Random_Sequence);Uxx=Uxx(Random_Sequence);Uxxx=Uxxx(Random_Sequence);

% Define the train data length
Train_Length=round(0.8*Data_Length);

% Define the training data
U_train=U(1:Train_Length,1);Ut_train=Ut(1:Train_Length,1);Utt_train=Utt(1:Train_Length,1);
Ux_train=Ux(1:Train_Length,1);Uxx_train=Uxx(1:Train_Length,1);Uxxx_train=Uxxx(1:Train_Length,1);

% Define the testing data
U_test=U(Train_Length+1:end,1);Ut_test=Ut(Train_Length+1:end,1);Utt_test=Utt(Train_Length+1:end,1);
Ux_test=Ux(Train_Length+1:end,1);Uxx_test=Uxx(Train_Length+1:end,1);Uxxx_test=Uxxx(Train_Length+1:end,1);

% Stack the test data for later use
Data_Test=[U_test,Ut_test,Utt_test,Ux_test,Uxx_test,Uxxx_test];

% Generate basic variables
syms u ut utt ux uxx uxxx f_x
Vars=[u ut utt ux uxx uxxx];

%% Now build your library to find the null space
% Original library
[SINDy_Data,SINDy_Struct]=SINDyLib_PDE_Ex1(U_train,Ut_train,Utt_train,Ux_train,Uxx_train,Uxxx_train);

% Use Null-Space method and set the tolerance
tol = 1e-15;

% Plotting option
pflag=0;

% Calculate the possible sparse null space vectors
[Xi, indTheta, lambdavec, numterms, errorv] = ADMpareto(SINDy_Data, tol, pflag);

%% Get the corresponding models
for j=1:size(Xi,2)
    % Calculate the correponding model, if the ut is not in the expression
    % then the symbolic solve will cause error, this suggests we can not
    % calculate the correct model.
    try
        digits 4
        PDE_RHS(j,1)=simplify(vpa(solve(cell2sym(SINDy_Struct)*Xi(:,j)==0,ut)));
        
        % Print the model we find
        fprintf('The %i possible model is: %s \n\n',j,char(PDE_RHS(j,1)))
        
        % Generate this function file
        Generate_PDE_RHS(PDE_RHS(j,1),Vars,strcat('TempFunctions/ParFor',num2str(j)));
        
        % Calculate the prediction error of the LHS
        L2_Error(j,1)=Get_Score(Ut_test,Data_Test,strcat('ParFor',num2str(j)));
        
        % Print score
        fprintf('\t Corresponding normalized l2 error is %i:\n',L2_Error(j,1))
    catch
        % Print the model we find
        fprintf('The %i possible model does not exist since ut is been ignored...\n\n',j)
        
        % Set corresponding error as NaN
        L2_Error(j,1)=NaN;
        
    end
end

%% Print the model that has the lowest l2 error
[minVal_Lambda,Index1]=min(L2_Error);

% Print the result
fprintf('\n\n\n The best model is:\n \t %s \n',PDE_RHS(Index1,1))

fprintf('\n The prediction error of this model is %i:\n',L2_Error(Index1,1))

%% Plot the LHS guess error
figure()
hold on
for jj=1:size(L2_Error,1)
    if isnan(L2_Error(jj,1))==0
        plot(jj,L2_Error(jj,1),'linewidth',2.5,'marker','o','color','green')
        set(gca,'FontSize',18);
        set(gcf,'Position',[100 100 600 400]);
        set(gcf,'PaperPositionMode','auto');
        grid on
    else
        plot(jj,0,'linewidth',2.5,'marker','x','color','black')
        set(gca,'FontSize',18);
        set(gcf,'Position',[100 100 600 400]);
        set(gcf,'PaperPositionMode','auto');
        grid on
    end
end












