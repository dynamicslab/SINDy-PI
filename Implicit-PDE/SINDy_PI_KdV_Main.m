%% This file is used to demonstrate that SINDy-PI could be used to discover the implicit PDE.
% Coded By: K
% Last Updated: 2019/06/18
%% Clear all
clc;clear;close all

%% Create folder and read function files
addpath('Functions')
[fld_status, fld_msg, fld_msgID]=mkdir('TempFunctions');
addpath('TempFunctions')

addpath('Functions')
addpath('TempFunctions')
addpath('./Datas')
Load_File='SimulationData_V_1.mat';
load(Load_File);

[fld_status, fld_msg, fld_msgID]=mkdir('SINDy_PI_Results');
addpath('SINDy_PI_Results')

% Define parameters
dx=x(2)-x(1);
PDE_Plot=1;

%% Now prepare the library data

% Reshape into a vector
U=real(reshape(u',[],1));
Ut=real(reshape(ut',[],1));
Ux=real(reshape(ux',[],1));
Uxx=real(reshape(uxx',[],1));
Uxxx=real(reshape(uxxx',[],1));
%% Plot the result
if PDE_Plot==1
    figure(1)
    hold on
    %
    surf(x,tspan,real(u))
    shading interp
    colormap(parula(100))
    %
    Spacer1=5;
    Spacer2=50;
    %
    ss=surf(x(1:Spacer1:size(x,2)),tspan(1:Spacer2:size(tspan,2)),real(u(1:Spacer2:size(tspan,2),1:Spacer1:size(x,2))))
    ss.EdgeColor='black'
    %
    set(gca,'FontSize',18);
    set(gcf,'Position',[100 100 600 400]);
    set(gcf,'PaperPositionMode','auto');
    grid on

    view(15,65)
end

%% Here, we prepare our data into testing set and training set
Data_Length=size(U,1);
Random_Sequence=randperm(Data_Length);

% Permute the original data
U=U(Random_Sequence);Ut=Ut(Random_Sequence);Utt=zeros(size(Ut));
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

%% Now build your LHS and RHS library
% Original library
[SINDy_Data,SINDy_Struct]=SINDyLib_DL_SINDy_mKdV(U_train,Ut_train,Utt_train,Ux_train,Uxx_train,Uxxx_train);

% LHS
[LHS_Data,LHS_Sym]=LHS_Guess_Implicit_PDE_mKdV(U_train,Ut_train,Utt_train,Ux_train,Uxx_train,Uxxx_train);
[LHS_Data_Test,~]=LHS_Guess_Implicit_PDE_mKdV(U_test,Ut_test,Utt_test,Ux_test,Uxx_test,Uxxx_test);

%% Now run SINDy-PI
% Define some parameter for sparese regression
N=20;
disp=0;
NormalizeLib=1;
dlambda=0.1;
sweep_len=100;

% Define some matrix to store the final result
Xi=cell(size(LHS_Sym,2),sweep_len);
PDE_RHS=cell(size(LHS_Sym,2),sweep_len);
L2_Error=zeros(size(LHS_Sym,2),sweep_len);

% We first loop through the LHS, then loop through the lambda
for iter=1:size(LHS_Sym,2)
    % Print the left hand side that we are testing
    fprintf('Testing the left hand side as %s:\n\n',char(LHS_Sym{iter}))
    
    % Get rid of LHS guess from the RHS library
    [RHS_Data,RHS_Sym]=ExcludeGuess_PDE(SINDy_Data,SINDy_Struct,LHS_Sym{iter});
    
    % Prepare variabels for parfor
    LHS_Dum=LHS_Sym{iter};
    LHS_Data_Dum=LHS_Data(:,iter);
    LHS_Data_Test_Dum=LHS_Data_Test(:,iter);
    
    for j=1:sweep_len
        % Print the left hand side that we are testing
        fprintf('\t Testing the lambda as %i:\n',j*dlambda)
        
        % SINDy-PI regression
        [Xi{iter,j},PDE_RHS{iter,j}] = sparsifyDynamics(RHS_Data,LHS_Data_Dum,LHS_Dum,j*dlambda,N,RHS_Sym,disp,NormalizeLib);
        
        % Solve for the PDE
        PDEs=solve(LHS_Dum==PDE_RHS{iter,j},ut);
        
        % Generate this function file
        Generate_PDE_RHS(PDEs,Vars,strcat('TempFunctions/ParFor',num2str(j)));
        
        % Calculate the prediction error of the LHS
        try
            L2_Error(iter,j)=Get_Score(Ut_test,Data_Test,strcat('ParFor',num2str(j)));
        catch
            % if the prediction error calculation throws an error it means
            % the symbolic sovling of the u_{t} throws empty equation
            L2_Error(iter,j)=NaN;
        end
        % Print score
        fprintf('\t Corresponding normalized l2 error is %i:\n',L2_Error(iter,j))
        
    end 
end

%% Print the model that has the lowest l2 error
for kk=1:size(L2_Error) 
    dummy=L2_Error(kk,:);
    [minVal_Lambda(kk,1),Index1(kk,1)]=min(dummy);
end

[minVal_Final,Index2]=min(minVal_Lambda);

% Print the result
digits(10)
fprintf('\n\n\n The best model is:\n \t %s \n',vpa(LHS_Sym{Index2}==PDE_RHS{Index2,Index1(Index2,1)}))

fprintf('\n The prediction error of this model is %i:\n',L2_Error(Index2,Index1(Index2,1)))

%% Store the best expression
Expression=vpa(LHS_Sym{Index2}==PDE_RHS{Index2,Index1(Index2,1)});

%% Save results
name=strcat('SINDy_PI_Results/',Load_File,'_Result.mat')
save(name,'Xi','PDE_RHS','sweep_len','dlambda','L2_Error','Expression','minVal_Final')
