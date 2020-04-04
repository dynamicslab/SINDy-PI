%% This file is used to test the final performance of the discovered SINDy-PI model.
% Please mannuly select and modify the initial condition as you wish.
% Coded By: K.Kahirman
% Last Updated: 2019/06/18
%% Clear all
clc;clear;close all

%% Create folder and read function files
addpath('Functions')
[fld_status, fld_msg, fld_msgID]=mkdir('TempFunctions');
addpath('TempFunctions')

%% We first generate the PDE data
% Define the time span
dt=0.01;T=20;
tspan=0:dt:T;

% Define the spatial domain
L=50;

% Define the discretization point
n=256;

% Whether you want to plot it or not
PDE_Plot=1;

% Whether compare?
Compare=0;

% Run simulation and gather data
[utsol1,usol1,x,k]=Example_1(dt,T,L,n,PDE_Plot,Compare);

[utsol2,usol2,x,k]=Example_1_Best_Model(dt,T,L,n,PDE_Plot);

%% Calculate Error
usol_error=abs(usol1-usol2)./abs(usol1);

%% Plot
close all
figure()
hold on

surf(x,tspan(1:1600),real(usol_error(:,1:1600)'),'FaceColor',[0.5 0.6 0.9],'EdgeColor','none','FaceLighting','none')
%
surf(x,tspan(1600:end),real(usol_error(:,1600:end)'),'FaceColor',[0.9 0.5 0.6],'EdgeColor','none','FaceLighting','none')
%
% surf(x,tspan,real(usol_error'))
% shading interp
% colormap(parula(100))
%
Spacer1=5;
Spacer2=50;
%
ss=surf(x(1:Spacer1:size(x,2)),tspan(1:Spacer2:size(tspan,2)),real(usol_error(1:Spacer1:size(x,2),1:Spacer2:size(tspan,2))'),'FaceAlpha',0)
ss.EdgeColor='black'
%
set(gca,'FontSize',24);
set(gcf,'Position',[100 100 600 400]);
set(gcf,'PaperPositionMode','auto');
grid on

%view(15,65)
view(42,35)
zlim([0,0.2])

figure()
hold on
%
surf(x,tspan(1:1600),real(usol_error(:,1:1600)'),'FaceColor',[0.5 0.6 0.9],'EdgeColor','none','FaceLighting','none')
%
surf(x,tspan(1600:end),real(usol_error(:,1600:end)'),'FaceColor',[0.9 0.5 0.6],'EdgeColor','none','FaceLighting','none')
%
%
% surf(x,tspan,real(usol_error'))
% shading interp
% colormap(parula(100))
%
Spacer1=5;
Spacer2=50;
%
ss=surf(x(1:Spacer1:size(x,2)),tspan(1:Spacer2:size(tspan,2)),real(usol_error(1:Spacer1:size(x,2),1:Spacer2:size(tspan,2))'),'FaceAlpha',0)
ss.EdgeColor='black'
%
set(gca,'FontSize',36);
set(gcf,'Position',[100 100 600 400]);
set(gcf,'PaperPositionMode','auto');
grid on

%view(15,65)
view(42,35)
%zlim([0,0.5])

% figure()
% pcolor(x,tspan,real(usol_error'))
% shading interp
% set(gca,'FontSize',18);
% set(gcf,'Position',[100 100 600 400]);
% set(gcf,'PaperPositionMode','auto');
% grid on
% colorbar










