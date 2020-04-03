%% This file is the main file of using SINDy-PI to discover the BZ-Reaction model
% Coded By: K
% Last Updated: 2019/06/25
%%
clc;close all;clear all;
%% Add folders to path
[fld_status, fld_msg, fld_msgID]=mkdir('TempFunctions');
[fld_status, fld_msg, fld_msgID]=mkdir('SINDy_PI_Results');
addpath('Datas')
addpath('Functions')
addpath('TempFunctions')
addpath('SINDy_PI_Results')
%% Now load the simulation file to get the simulation data
fprintf('\n Loading the file...')
tic
Load_File='Simulation_BZ_0.1_Seconds_128_Grid.mat';
load(strcat('Datas\',Load_File))
fprintf('\n Loading finished, used %d seconds!\n',toc)
%% Now build the base library
fprintf('\n Reshaping the vector...')
% Reshape into a vector
tic
R=real(reshape(r,[],1));Rt=real(reshape(r_t,[],1));Rx=real(reshape(r_x,[],1));Rxx=real(reshape(r_xx,[],1));...
    Ry=real(reshape(r_y,[],1));Ryy=real(reshape(r_yy,[],1));
Z=real(reshape(z,[],1));Zt=real(reshape(z_t,[],1));Zx=real(reshape(z_x,[],1));Zxx=real(reshape(z_xx,[],1));...
    Zy=real(reshape(z_y,[],1));Zyy=real(reshape(z_yy,[],1));
S=real(reshape(s,[],1));St=real(reshape(s_t,[],1));Sx=real(reshape(s_x,[],1));Sxx=real(reshape(s_xx,[],1));...
    Sy=real(reshape(s_y,[],1));Syy=real(reshape(s_yy,[],1));
U=real(reshape(u,[],1));Ut=real(reshape(u_t,[],1));Ux=real(reshape(u_x,[],1));Uxx=real(reshape(u_xx,[],1));...
    Uy=real(reshape(u_y,[],1));Uyy=real(reshape(u_yy,[],1));
%
fprintf('\n Reshape finished, used %d seconds!\n',toc)

%% Here, we prepare our data into testing set and training set
fprintf('\n Creating the training and testing data...')

tic
Data_Length=size(U,1);
Random_Sequence=randperm(Data_Length);

% Permute the original data
R=R(Random_Sequence);Rt=Rt(Random_Sequence);Rx=Rx(Random_Sequence);Rxx=Rxx(Random_Sequence);Ry=Ry(Random_Sequence);Ryy=Ryy(Random_Sequence);
Z=Z(Random_Sequence);Zt=Zt(Random_Sequence);Zx=Zx(Random_Sequence);Zxx=Zxx(Random_Sequence);Zy=Zy(Random_Sequence);Zyy=Zyy(Random_Sequence);
S=S(Random_Sequence);St=St(Random_Sequence);Sx=Sx(Random_Sequence);Sxx=Sxx(Random_Sequence);Sy=Sy(Random_Sequence);Syy=Syy(Random_Sequence);
U=U(Random_Sequence);Ut=Ut(Random_Sequence);Ux=Ux(Random_Sequence);Uxx=Uxx(Random_Sequence);Uy=Uy(Random_Sequence);Uyy=Uyy(Random_Sequence);

% Define the train data length
Train_Length=round(0.8*Data_Length);

% Define the training data
R_train=R(1:Train_Length,1);Rt_train=Rt(1:Train_Length,1);
Rx_train=Rx(1:Train_Length,1);Rxx_train=Rxx(1:Train_Length,1);
Ry_train=Ry(1:Train_Length,1);Ryy_train=Ryy(1:Train_Length,1);
%
Z_train=Z(1:Train_Length,1);Zt_train=Zt(1:Train_Length,1);
Zx_train=Zx(1:Train_Length,1);Zxx_train=Zxx(1:Train_Length,1);
Zy_train=Zy(1:Train_Length,1);Zyy_train=Zyy(1:Train_Length,1);
%
S_train=S(1:Train_Length,1);St_train=St(1:Train_Length,1);
Sx_train=Sx(1:Train_Length,1);Sxx_train=Sxx(1:Train_Length,1);
Sy_train=Sy(1:Train_Length,1);Syy_train=Syy(1:Train_Length,1);
%
U_train=U(1:Train_Length,1);Ut_train=Ut(1:Train_Length,1);
Ux_train=Ux(1:Train_Length,1);Uxx_train=Uxx(1:Train_Length,1);
Uy_train=Uy(1:Train_Length,1);Uyy_train=Uyy(1:Train_Length,1);

% Define the testing data
R_test=R(Train_Length+1:end,1);Rt_test=Rt(Train_Length+1:end,1);
Rx_test=Rx(Train_Length+1:end,1);Rxx_test=Rxx(Train_Length+1:end,1);
Ry_test=Ry(Train_Length+1:end,1);Ryy_test=Ryy(Train_Length+1:end,1);
%
Z_test=Z(Train_Length+1:end,1);Zt_test=Zt(Train_Length+1:end,1);
Zx_test=Zx(Train_Length+1:end,1);Zxx_test=Zxx(Train_Length+1:end,1);
Zy_test=Zy(Train_Length+1:end,1);Zyy_test=Zyy(Train_Length+1:end,1);
%
S_test=S(Train_Length+1:end,1);St_test=St(Train_Length+1:end,1);
Sx_test=Sx(Train_Length+1:end,1);Sxx_test=Sxx(Train_Length+1:end,1);
Sy_test=Sy(Train_Length+1:end,1);Syy_test=Syy(Train_Length+1:end,1);
%
U_test=U(Train_Length+1:end,1);Ut_test=Ut(Train_Length+1:end,1);
Ux_test=Ux(Train_Length+1:end,1);Uxx_test=Uxx(Train_Length+1:end,1);
Uy_test=Uy(Train_Length+1:end,1);Uyy_test=Uyy(Train_Length+1:end,1);

% Stack the training data to generate a library
Data_Train=[R_train Rt_train Rx_train Rxx_train Ry_train Ryy_train Z_train Zt_train Zx_train Zxx_train Zy_train Zyy_train...
    S_train St_train Sx_train Sxx_train Sy_train Syy_train U_train Ut_train Ux_train Uxx_train Uy_train Uyy_train];

% Stack the test data for later use
Data_Test=[R_test Rt_test Rx_test Rxx_test Ry_test Ryy_test Z_test Zt_test Zx_test Zxx_test Zy_test Zyy_test...
    S_test St_test Sx_test Sxx_test Sy_test Syy_test U_test Ut_test Ux_test Uxx_test Uy_test Uyy_test];

% Generate basic variables
syms r rt rx rxx ry ryy z zt zx zxx zy zyy ...
    s st sx sxx sy syy u ut ux uxx uy uyy
Vars=[r rt rx rxx ry ryy z zt zx zxx zy zyy s st sx sxx sy syy u ut ux uxx uy uyy];
%
fprintf('\n Training and testing data creation finished, used %d seconds!\n',toc)


%% Now run SINDy-PI
fprintf('\n Using SINDy-PI on the data...\n\n')

fprintf('\t \n Setting up the parameters...\n\n')
% Define some parameter for sparese regression
N=20;
disp=0;
NormalizeLib=1;
dlambda=0.1;
sweep_len=10;
Eqn_Num=4; % Define how many equations you have

% Define some matrix to store the final result
Xi=cell(Eqn_Num,size(3,2),sweep_len);
PDE_RHS=cell(Eqn_Num,size(3,2),sweep_len);
L2_Error=zeros(Eqn_Num,size(3,2),sweep_len);

tic
for WhichEqn=1:Eqn_Num
    fprintf('\t \n Calculating the %i equation...\n\n',WhichEqn)
    
    % Original library
    [SINDy_Data,SINDy_Struct]=SINDyLib_DL_SINDy_BZ_Reaction(Data_Train,Vars,WhichEqn);
    
    % LHS
    [LHS_Data,LHS_Sym]=SINDyLib_DL_SINDy_BZ_Reaction(Data_Train,Vars,WhichEqn);
    [LHS_Data_Test,~]=SINDyLib_DL_SINDy_BZ_Reaction(Data_Test,Vars,WhichEqn);
    
    % We first loop through the LHS, then loop through the lambda
    for iter=1:size(LHS_Sym,2)
        % Print the left hand side that we are testing
        fprintf('\t\t Testing the left hand side as %s:\n\n',char(LHS_Sym{iter}))
        
        % Get rid of LHS guess from the RHS library
        [RHS_Data,RHS_Sym]=ExcludeGuess_PDE(SINDy_Data,SINDy_Struct,LHS_Sym{iter});
        
        % Prepare variabels for parfor
        LHS_Dum=LHS_Sym{iter};
        LHS_Data_Dum=LHS_Data(:,iter);
        LHS_Data_Test_Dum=LHS_Data_Test(:,iter);
        
        %parfor j=1:sweep_len
        for j=1:sweep_len
            % Print the left hand side that we are testing
            fprintf('\t\t Testing the lambda as %i...\n',j*dlambda)
            
            % SINDy-PI regression
            [Xi{WhichEqn,iter,j},PDE_RHS{WhichEqn,iter,j}] = sparsifyDynamics(RHS_Data,LHS_Data_Dum,LHS_Dum,j*dlambda,N,RHS_Sym,disp,NormalizeLib);
            %lambda=0.1;
            %[Xi{WhichEqn,iter,j},PDE_RHS{WhichEqn,iter,j}] = sparsifyDynamics(RHS_Data,LHS_Data_Dum,LHS_Dum,lambda,N,RHS_Sym,disp,NormalizeLib);
            
            % Generate this function file
            Generate_PDE_RHS(PDE_RHS{WhichEqn,iter,j},Vars,strcat('TempFunctions/ParFor',num2str(j)));
            
            % Calculate the prediction error of the LHS
            L2_Error(WhichEqn,iter,j)=Get_Score(LHS_Data_Test_Dum,Data_Test,strcat('ParFor',num2str(j)));
            
            % Print score
            fprintf('\t\t Corresponding normalized l2 error is %i:\n',L2_Error(WhichEqn,iter,j))
            
        end
    end
    
    % Now output the best model of current equation, print the model that has the lowest l2 error
    
    for kk=1:size(LHS_Sym,2)
        dummy=L2_Error(WhichEqn,kk,:);
        [minVal_Lambda(WhichEqn,kk),Index1(WhichEqn,kk)]=min(dummy);
    end
    
    [minVal_Final(WhichEqn,1),Index2(WhichEqn,1)]=min(minVal_Lambda(WhichEqn,:));
    
    % Print the result
    digits(4)
    fprintf('\n\n\n  \t\t\t  The best model is: %s \n \t  \n',...
        vpa(LHS_Sym{Index2(WhichEqn,1)}==PDE_RHS{WhichEqn,Index2(WhichEqn,1),Index1(WhichEqn,Index2(WhichEqn,1))}))
    
    fprintf('\n \t\t\t The prediction error of this model is: %i\n',...
        L2_Error(WhichEqn,Index2(WhichEqn,1),Index1(WhichEqn,Index2(WhichEqn,1))))
    
    % Store the best expression
    Expression(WhichEqn,1)=vpa(LHS_Sym{Index2(WhichEqn,1)}==PDE_RHS{WhichEqn,Index2(WhichEqn,1),Index1(WhichEqn,Index2(WhichEqn,1))});

end


%% Save results
Save_Final_Result=1;
if Save_Final_Result==1
    fprintf('\n Saving the result to the folder...')
    tic
    name=strcat('SINDy_PI_Results/',Load_File,'_Result.mat')
    save(name,'Xi','PDE_RHS','sweep_len','dlambda','L2_Error','Expression','minVal_Final')
    fprintf('\n Saving process finished, used %d seconds!\n',toc)
end


















