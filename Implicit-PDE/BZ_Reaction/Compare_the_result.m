%% This file will be used to compare identified PDE of the SINDy-PI and ground truth.
% Code By: K
% Last Updated: 2019/06/24
%% Clear all
close all;clear all;clc;
%% Add path
[fld_status, fld_msg, fld_msgID]=mkdir('Datas');
[fld_status, fld_msg, fld_msgID]=mkdir('Figures');
addpath('Datas')
addpath('Figures')
addpath('Functions')
%% Define parameters for the simulation
% Define diffusion rate
Dx=0.01;Dz=0.01;Ds=1;Du=1;

% Others
q=0.1;f=1.5;ksi=0.3;alpha=0.3;beta=0.26;gama=0.4;ksi2=1.5;ksi3=0.003;phi=0;

% Define the time horizon
dt=0.001;T=1;
tspan=0:dt:T;

% Define the spatial domain
L=20; % Total length of each dimension
n=128; % Discretization point of each dimension
N=n*n; % Total points used, n^2

x2=linspace(-L/2,L/2,n+1); x=x2(1:n);
y=x;dx=x(2)-x(1);

kx=(2*pi/L)*[0:(n/2-1) -n/2:-1];
ky=kx;

% Get n-dimensional grid
[X,Y]=meshgrid(x,y);
[KX,KY]=meshgrid(kx,ky);
K2=KX.^2+KY.^2; K22=reshape(K2,N,1);
Kx=reshape(KX,N,1);Ky=reshape(KY,N,1);Kxx=reshape(KX,N,1);Kyy=reshape(KY,N,1);

% Create a time matrix
r_t=zeros(n,n,length(tspan));r=zeros(n,n,length(tspan));r_x=zeros(n,n,length(tspan));r_xx=zeros(n,n,length(tspan));r_y=zeros(n,n,length(tspan));r_yy=zeros(n,n,length(tspan));
z_t=zeros(n,n,length(tspan));z=zeros(n,n,length(tspan));z_x=zeros(n,n,length(tspan));z_xx=zeros(n,n,length(tspan));z_y=zeros(n,n,length(tspan));z_yy=zeros(n,n,length(tspan));
s_t=zeros(n,n,length(tspan));s=zeros(n,n,length(tspan));s_x=zeros(n,n,length(tspan));s_xx=zeros(n,n,length(tspan));s_y=zeros(n,n,length(tspan));s_yy=zeros(n,n,length(tspan));
u_t=zeros(n,n,length(tspan));u=zeros(n,n,length(tspan));u_x=zeros(n,n,length(tspan));u_xx=zeros(n,n,length(tspan));u_y=zeros(n,n,length(tspan));u_yy=zeros(n,n,length(tspan));
r_lap=zeros(n,n,length(tspan));z_lap=zeros(n,n,length(tspan));s_lap=zeros(n,n,length(tspan));u_lap=zeros(n,n,length(tspan));
%% Set a Guassian for the initial condition
r0=cos(sqrt(X.^2+Y.^2)).^2+GaussianFilter(X,Y,0.05,0.05,0,0)+GaussianFilter(X,Y,1,0.01,0,0)+GaussianFilter(X,Y,0.01,1,0,0)+GaussianFilter(X,Y,0.25,0.25,5,5)+GaussianFilter(X,Y,0.25,0.25,5,-5)+GaussianFilter(X,Y,0.25,0.25,-5,5)+GaussianFilter(X,Y,0.25,0.25,-5,-5)+0.1;
z0=sin(0.1*sqrt(X.^2+Y.^2)).^2+GaussianFilter(X,Y,0.05,0.01,0,0)+GaussianFilter(X,Y,0.25,0.25,5,5)+GaussianFilter(X,Y,0.25,0.25,5,-5)+GaussianFilter(X,Y,0.25,0.25,-5,5)+GaussianFilter(X,Y,0.25,0.25,-5,-5)+0.1;
s0=sin(sqrt(X.^2+Y.^2)).^2+GaussianFilter(X,Y,0.5,0.5,0,5)+GaussianFilter(X,Y,0.5,0.5,5,0)+GaussianFilter(X,Y,0.5,0.5,0,-5)+GaussianFilter(X,Y,0.5,0.5,-5,0)+0.1;
u0=cos(sqrt(X.^2+Y.^2)).^2+GaussianFilter(X,Y,0.5,0.01,0,0)+GaussianFilter(X,Y,0.01,0.5,0,0)+GaussianFilter(X,Y,0.25,0.25,5,5)+GaussianFilter(X,Y,0.25,0.25,5,-5)+GaussianFilter(X,Y,0.25,0.25,-5,5)+GaussianFilter(X,Y,0.25,0.25,-5,-5)+0.01;

%% Plot the initial conditions
figure(1)
surf(x,y,r0)
set(gca,'FontSize',18);
set(gcf,'Position',[100 100 600 600]);
set(gcf,'PaperPositionMode','auto');
view(-18,35)
box('off')
axis('off')

figure(2)
surf(x,y,z0)
set(gca,'FontSize',18);
set(gcf,'Position',[100 100 600 600]);
set(gcf,'PaperPositionMode','auto');
view(-18,35)
box('off')
axis('off')

figure(3)
surf(x,y,s0)
set(gca,'FontSize',18);
set(gcf,'Position',[100 100 600 600]);
set(gcf,'PaperPositionMode','auto');
view(-18,35)
box('off')
axis('off')

figure(4)
surf(x,y,u0)
set(gca,'FontSize',18);
set(gcf,'Position',[100 100 600 600]);
set(gcf,'PaperPositionMode','auto');
view(-18,35)
box('off')
axis('off')
%%
close all
% Transfer the initial guess to fourier domain
xzsu0t=[reshape(fft2(r0),1,N) reshape(fft2(z0),1,N) reshape(fft2(s0),1,N) reshape(fft2(u0),1,N)].';

%% Simulate the original equation
tic
opts = odeset('RelTol',1e-12,'AbsTol',1e-13);
NeedDev=0;
[time,xzsusol]=ode45(@(time,xzsut)BZ_Reaction_PDE(time,xzsut,Kx,Kxx,Ky,Kyy,K22,n,N,Dx,Dz,Ds,Du,q,f,ksi,alpha,beta,gama,ksi2,ksi3,phi,NeedDev),tspan,xzsu0t,opts);
toc

% Simulate the DL-SINDy discovered equation
tic
opts = odeset('RelTol',1e-12,'AbsTol',1e-13);
NeedDev=0;
[time,xzsusol_DL]=ode45(@(time,xzsut)BZ_Reaction_DL_SINDy_PDE(time,xzsut,Kx,Kxx,Ky,Kyy,K22,n,N,NeedDev),tspan,xzsu0t,opts);
toc

%% Extract values
% After simulation, we extract the derivative
NeedDev=1;

for pin=1:length(tspan)
    % Get the derivative data
    [~,r_t(:,:,pin),z_t(:,:,pin),s_t(:,:,pin),u_t(:,:,pin),r(:,:,pin),z(:,:,pin),s(:,:,pin),u(:,:,pin),r_x(:,:,pin),z_x(:,:,pin),s_x(:,:,pin),u_x(:,:,pin),r_y(:,:,pin),z_y(:,:,pin),s_y(:,:,pin),u_y(:,:,pin),r_xx(:,:,pin),z_xx(:,:,pin),s_xx(:,:,pin),u_xx(:,:,pin),r_yy(:,:,pin),z_yy(:,:,pin),s_yy(:,:,pin),u_yy(:,:,pin),r_lap(:,:,pin),z_lap(:,:,pin),s_lap(:,:,pin),u_lap(:,:,pin)]=...
        BZ_Reaction_PDE(0,xzsusol(pin,:).',Kx,Kxx,Ky,Kyy,K22,n,N,Dx,Dz,Ds,Du,q,f,ksi,alpha,beta,gama,ksi2,ksi3,phi,NeedDev);

    % Do the same for the DL-SINDy
    [~,~,~,~,~,r_DL(:,:,pin),z_DL(:,:,pin),s_DL(:,:,pin),u_DL(:,:,pin),~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~]=...
        BZ_Reaction_SINDy_PI_PDE(0,xzsusol_DL(pin,:).',Kx,Kxx,Ky,Kyy,K22,n,N,NeedDev);
end

%% Plot the result for original data
for pin=2:4
    if pin==1
        Animate=r;
        limits=[0 1];
    elseif pin==2
        Animate=z;
        limits=[0 1];
    elseif pin==3
        Animate=s;
        limits=[0 1];
    elseif pin==4
        Animate=u;
        limits=[0 1];
    end
    % Normalize it
    Dummy=Animate(:,:,1);
    Max=max(max(Dummy))
    Animate=Animate/Max;
    %
    for j=1:200:length(tspan)
        clf
        surf(x,y,Animate(:,:,j))
        colormap(cool)
        %shading interp
        %caxis([0,1])
        %Make plot propotional
        %zlim([0,1.1])
        %set(gca,'FontSize',18);
        set(gcf,'Position',[100 100 600 600]);
        set(gcf,'PaperPositionMode','auto');
        view(-18,35)
        box('off')
        axis('off')
        drawnow limitrate
    end
end

%% Plot the result for DL-SINDy simulation
for pin=2:4
    if pin==1
        Animate=r_DL;
        limits=[0 1];
    elseif pin==2
        Animate=z_DL;
        limits=[0 1];
    elseif pin==3
        Animate=s_DL;
        limits=[0 1];
    elseif pin==4
        Animate=u_DL;
        limits=[0 1];
    end
    % Normalize it
    Dummy=Animate(:,:,1);
    Max=max(max(Dummy))
    Animate=Animate/Max;
    %
    for j=1:200:length(tspan)
        clf
        surf(x,y,Animate(:,:,j))
        colormap(cool)
        %caxis([0,1])
        %Make plot propotional
        %zlim([0,1.1])
        %set(gca,'FontSize',18);
        set(gcf,'Position',[100 100 600 600]);
        set(gcf,'PaperPositionMode','auto');
        view(-18,35)
        box('off')
        axis('off')
        drawnow limitrate
    end
end










