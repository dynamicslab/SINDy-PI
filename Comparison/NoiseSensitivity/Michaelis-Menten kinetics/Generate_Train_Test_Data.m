%% This file will be used to generate the training and testing data. This data
% will then be provided to the SINDy-PI and implicit SINDy to calculate the
% system equation.
% Coded By: K
% Last Updated: 2019/07/11
%% Clear
clc;clear all;close all;
%% First we generate a folder to store all the simulation data
[status,msg] = mkdir('Datas');
addpath('Datas')
addpath('Functions')

%% Now generate data
% Define simulation parameter
dt=0.1; T=5;

% Define how many train and test data you need, the ratio we use here is
% 8:2
init_num_Train=2400;
Gap_Train=200;
init_num_Test=600;
Gap_Test=50;


% Determine how many iteration you need
Final_pin=30;

% How many noise level
percent_start=1;d_percent=1;percent_end=24;

% How many noisy data generation method
Num_Method=3;

for iter=1:Final_pin
    fprintf('Now calculating iteration %d...\n\n',iter)
    percent=0;
    
    for percent_iter=percent_start:d_percent:percent_end
        
        % Set the pin one step forward
        percent=percent+1;
        
        % Set the noise level
        if percent_iter==1
            noise_level(percent,1)=0;
        elseif percent_iter==2
            noise_level(percent,1)=1e-7;
        elseif percent_iter==3
            noise_level(percent,1)=5e-7;
        elseif percent_iter==4
            noise_level(percent,1)=1e-6;
        elseif percent_iter==5
            noise_level(percent,1)=5e-6;
        elseif percent_iter==6
            noise_level(percent,1)=1e-5;
        elseif percent_iter==7
            noise_level(percent,1)=5e-5;
        elseif percent_iter==8
            noise_level(percent,1)=1e-4;
        elseif percent_iter==9
            noise_level(percent,1)=5e-4;
        elseif percent_iter==10
            noise_level(percent,1)=1e-3;
        elseif percent_iter<=14
            noise_level(percent,1)=2e-3*(percent_iter-10);
        elseif percent_iter==15
            noise_level(percent,1)=1e-2;
        elseif percent_iter<=19
            noise_level(percent,1)=2e-2*(percent_iter-15);
        elseif percent_iter==20
            noise_level(percent,1)=1e-1;
        elseif percent_iter<=24
            noise_level(percent,1)=1e-1+1e-1*(percent_iter-20);
        end
        
        fprintf('\t Calculating noise level as %d...\n',noise_level(percent,1))
        
        for Denoise_Method=1:Num_Method
            % Generate data
            [Data,dData,iter_length]=Generate_Noisy_Data(noise_level(percent,1),init_num_Train,dt,T,Denoise_Method,Gap_Train);
            [xt_test,dxt_test,iter_length]=Generate_Noisy_Data(noise_level(percent,1),init_num_Test,dt,T,Denoise_Method,Gap_Test);
            % Save it
            File_Name=strcat('Datas/Data_Iter_',num2str(iter),'_NoiseLevel_',num2str(noise_level(percent,1)),'_NoiseMethod_',num2str(Denoise_Method),'_GoodLuck.mat');
            save(File_Name,'Data','dData','xt_test','dxt_test')  
        end
    end
    
end



