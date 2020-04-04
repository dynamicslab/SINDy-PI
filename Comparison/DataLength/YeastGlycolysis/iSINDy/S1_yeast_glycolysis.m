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
% Find equations for yeast glycolysis state variable 1
%% 
clc;clear all;close all;

addpath('./utils');
addpath('./bioutils');

% define libarary parameters
laurentorder = 0;
polyorder = 6;
usesine = 0;
dyorder = 1;

% search for equation for the 1st state-variable
% We do the ADM method for only one value of lambda 
% see 2nd state varaiable for full sweep.
% clear results for other state variables
clear Theta Thetastring Xi indTheta lambdavec numterms errorv indopt
clear indTheta1 Xi1 numterms1 nT

% Here we load the previously simulated data
load('TrainingData.mat')

% Define the new data length
percent=0.1;new_length=round(percent*length(xt));

% Shuffel the original data
Sequence=randperm(size(xt,1));
xt=xt(Sequence,:);
dxt=dxt(Sequence,:);

% Assign the value to the new variables
Data=xt(1:new_length,:);dData=dxt(1:new_length,:);

% Define the number of states
n=size(Data,2);

% pool Data  (i.e., build library of nonlinear time series)
[Theta, Thetastring] = poolDatady(Data,n,polyorder,usesine, laurentorder, dData(:,1), dyorder);
% %initial lambda value, which is the value used for soft thresholding in ADM

%%
tol = 2e-3;

lambda = 5e-3;

jj = 1; % counter
num= 1; % initialize the number of nonzero terms found for the lambda
errorvec= 0;

MaxIter = 1e2;

% for now calculate null space using null function
nT = null(Theta);
plottag=2;
[indTheta1, Xi1, numterms1] = ADMinitvary(nT,lambda,MaxIter,tol, plottag);

Thetastring(indTheta1)'
n0Xi = Xi1(Xi1~=0); % terms need to be rearranged to recover coefficients
n0Xi/n0Xi(end)
%%
save('Results/1st_state_variable.mat')