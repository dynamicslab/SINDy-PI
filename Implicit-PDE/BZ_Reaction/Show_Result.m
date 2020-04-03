%% This file is used to read and plot the final result of the SINDy-PI
% Coded By: K
% Last Updated: 2019/06/26
%%
clc;clear all;close all;
% Add folders
addpath('SINDy_PI_Results')

% load file
load('Simulation_BZ_0.1_Seconds_128_Grid.mat_Result.mat')

%% Show the final PDE
digits(5)
vpa(Expression)

syms rt zt st ut;
State={rt,zt,st,ut};
for i=1:size(Expression,1)
    Eq(i,1)=simplify(solve(Expression(i),State{1,i}));
end

Eq




