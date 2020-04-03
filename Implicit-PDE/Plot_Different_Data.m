%% This file will be used to plot the different simulation data
% Coded By: K
% Last Updated: 2019/06/23
%%
clc;clear all; close all;

%% Add path
addpath('./Datas')

%% Load data
load('SimulationData_V_1.mat')

%% Plot
close all
% figure(1)
% waterfall(x,tspan,real(u))
% colormap([0 0 0]);
% view(42,55)
% set(gca,'FontSize',18);
% set(gcf,'Position',[100 100 600 400]);
% set(gcf,'PaperPositionMode','auto');
% grid on

figure(2)
Spacer1=10;
Spacer2=15;
h1=surf(x,tspan,real(u))
set(h1,'LineStyle','none')
view(42,55)
set(gca,'FontSize',40);
set(gcf,'Position',[100 100 600 400]);
set(gcf,'PaperPositionMode','auto');
grid on

hold on
h2=surf(x(1:Spacer1:end),tspan(1:Spacer2:end),real(u(1:Spacer2:length(tspan),1:Spacer1:length(x))))
set(h2,'LineStyle','-')

% figure(3)
% pcolor(x,tspan,real(u))
% shading interp
% set(gca,'FontSize',18);
% set(gcf,'Position',[100 100 600 400]);
% set(gcf,'PaperPositionMode','auto');
% grid on




