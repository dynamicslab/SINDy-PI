%% This file will be used to plot the different result of SINDy-PI
% Coded By: K
% Last Updated: 2019/06/23
%%
clc;clear all; close all;

%% Add path
addpath('.\SINDy_PI_Results')

%% Load data
Files=dir('.\SINDy_PI_Results\*.mat');

Error=zeros(length(Files),1);
for k=1:length(Files)
   FileNames=Files(k).name;
   load(FileNames)
   Error(k,1)=minVal_Final;
   digits(4)
   Expression
end


%% Plot
close all
g0=[0;0.0001;0.001;0.01;0.1;1];
figure(1)
loglog(g0,Error,'marker','o','linewidth',2.5)
set(gca,'FontSize',18);
set(gcf,'Position',[100 100 600 400]);
set(gcf,'PaperPositionMode','auto');
grid on


