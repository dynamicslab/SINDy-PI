%% This file will be used to plot the performance comparison of SINDy-PI and PDE-Find
% Coded By: K
% Last Updated: 2019/06/23
%%
clc;clear all; close all;

%% Add path
addpath('./SINDy_PI_Results')
addpath('./PDE_Find_Results')
%% Load data
Files=dir('./SINDy_PI_Results/*.mat');

Error_DL=zeros(length(Files),1);
for k=1:length(Files)
   FileNames=strcat('SINDy_PI_Results/',Files(k).name);
   load(FileNames)
   Error_DL(k,1)=minVal_Final;
   digits(4)
   Expression
end

%%
Files=dir('./PDE_Find_Results/*.mat');

Error_PDE_Find=zeros(length(Files),1);
for k=1:length(Files)
   FileNames=strcat('PDE_Find_Results/',Files(k).name);
   load(FileNames)
   Error_PDE_Find(k,1)=minVal_Final;
   digits(4)
   Expression
end

%% Plot
close all
% Note that the first point should have zero noise, we set it up as 1e-5
% for the convenience of plotting
g0=[1e-5;0.0001;0.001;0.01;0.1;1];
figure(1)
hold on
set(gca,'FontSize',18);
set(gcf,'Position',[100 100 800 150]);
set(gcf,'PaperPositionMode','auto');
grid on

plot(g0,Error_PDE_Find,'linewidth',3,'linestyle','-','color','red')
scatter(g0,Error_PDE_Find,100,'filled','MarkerEdgeColor',[0 0 0],'MarkerFaceColor',1*[1 0 0])

plot(g0,Error_DL,'linewidth',3,'linestyle','-','color','blue')
scatter(g0,Error_DL,100,'filled','MarkerEdgeColor',[0 0 0],'MarkerFaceColor',1*[0 0 1])

xticks([1e-5 1e-4 1e-3 1e-2 1e-1 1])
xticklabels({'0','10^{-4}','10^{-3}','10^{-2}','10^{-1}','1'})

box('on')
set(gca, 'XScale', 'log', 'YScale', 'log');



