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
% Find equations for yeast glycolysis state variable 5
%% 
clc;clear all;close all;

addpath('./utils');
addpath('./bioutils');

% define libarary parameters
laurentorder = 0;
polyorder = 3;
usesine = 0;
dyorder = 1;

% clear variables from other solutions.
clear Theta Thetastring Xi indTheta lambdavec numterms errorv indopt
clear indTheta1 Xi1 numterms1 nT

% Here we load the previously simulated data
load('TrainingData.mat')

% Define the new data length
percent=0.003;new_length=round(percent*length(xt));

% Shuffel the original data
Sequence=randperm(size(xt,1));
xt=xt(Sequence,:);
dxt=dxt(Sequence,:);

% Assign the value to the new variables
Data=xt(1:new_length,:);dData=dxt(1:new_length,:);

% Define the number of states
n=size(Data,2);

% pool Data  (i.e., build library of nonlinear time series)
[Theta, Thetastring] = poolDatady(Data,n,polyorder,usesine, laurentorder, dData(:,5), dyorder);
% %initial lambda value, which is the value used for soft thresholding in ADM

tol = 2e-3;
lambda = 8e-3;

jj = 1; % counter
num= 1; % initialize the number of nonzero terms found for the lambda
errorvec= 0;

MaxIter = 1e3;

% for now calculate null space using null function
nT = null(Theta);

% Define the parameters for plooting 
plottag=2;

[indTheta1, Xi1, numterms1] = ADMinitvary(nT,lambda,MaxIter,tol, plottag);
Thetastring(indTheta1)'
n0Xi = Xi1(Xi1~=0);
n0Xi/n0Xi(end)
save('Results/5th_state_variable.mat')
