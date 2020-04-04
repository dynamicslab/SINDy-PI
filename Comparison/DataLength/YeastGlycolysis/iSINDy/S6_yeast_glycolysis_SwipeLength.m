% Originally Coded By: Niall Mangan
% Copyright 2017, All Rights Reserved
% Code by Niall Mangan for paper "Inferring biological networks by sparse
% identification of nonlinear dynamics"
% by N. M. Mangan S. L. Brunton, J. L. Proctor, and J. N. Kutz
% Original Code: https://github.com/niallmm/iSINDy
%%
% This file is modified to compare the performance of i-SINDy on
% insufficient data.
% Modified by:K, 2019/07/16
% Find equations for yeast glycolysis state variable 6
%%
clc;clear all;close all;

addpath('./utils');
addpath('./bioutils');
FolderName='Results2';

% define libarary parameters
laurentorder = 0;
polyorder = 6;
usesine = 0;
dyorder = 1;

% clear variables from other solutions.
clear Theta Thetastring Xi indTheta lambdavec numterms errorv indopt
clear indTheta1 Xi1 numterms1 nT

% Here we load the previously simulated data
load('TrainingData.mat')

%% Define some parameters for the implicit SINDy
% Initial lambda value, which is the value used for soft thresholding in ADM
% To make implicit SINDy work, you need to carefully tune these two
% parameters. The below values are the best one that I found.
tol = 1e-4;
lambda = 3e-4;

% counter
jj = 1;

% initialize the number of nonzero terms found for the lambda
num= 1;
errorvec= 0;

% Define maximum iteration you need
MaxIter = 1000;

% Plot option, we set it to zero
plottag=0;

% Define the for loop
percent_start=0.3;
percent_end=1;
d_percent=0.1;

% Define how many results you want to get for each percentage
N_Iter=20;

for iter=1:N_Iter
    for percent=percent_start:d_percent:percent_end
        fprintf('\n\n\t Using %i percent of the data...\n',percent*100)
        
        fprintf('\n\t\t Calculating for the %i time...\n',iter)
        % Get the data length
        new_length=round(percent*length(xt));
        
        % Shuffel the original data
        Sequence=randperm(size(xt,1));
        xt_dum=xt(Sequence,:);
        dxt_dum=dxt(Sequence,:);
        
        % Assign the value to the new variables
        Data=xt_dum(1:new_length,:);dData=dxt_dum(1:new_length,:);
        
        % Define the number of states
        n=size(Data,2);
        
        % pool Data  (i.e., build library of nonlinear time series)
        [Theta, Thetastring] = poolDatady(Data,n,polyorder,usesine, laurentorder, dData(:,6), dyorder);
        
        % for now calculate null space using null function
        nT = null(Theta);
        
        % Get the sparse vector
        tic
        [indTheta1, Xi1, numterms1] = ADMinitvary(nT,lambda,MaxIter,tol, plottag);
        fprintf('\t\t\t Calculation finished! Used %d seconds.\n',toc)
        
        % Save the calculation result of current iteration
        fprintf('\n\t\t\t\t Saving the result...\n')
        tic
        cc=clock;
        ResultName=strcat(FolderName,'/implicit_SINDY_Data_Length_',num2str(percent*100),'_','Iter',num2str(iter),'_',num2str(cc(3)),'_',num2str(cc(4)),'_',num2str(cc(5)),'_',num2str(round(cc(6))),'.mat');
        save(ResultName,'Xi1','indTheta1','percent')
        fprintf('\n\t\t\t\t Saving finished! Using %i seconds...\n',toc)
        
    end
end
