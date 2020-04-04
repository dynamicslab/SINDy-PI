%% This file is used to demonstrate that SINDy-PI could be used to discover the implicit PDE.
% Coded By: K
% Last Updated: 2019/06/18
%% Clear all
clc;clear;close all
%% We first generate the PDE data
% Define the time span
dt=0.01;T=15;
tspan=0:dt:T;

% Define the spatial domain
L=50;

% Define the discretization point
n=256;

% Whether you want to plot it or not
PDE_Plot=0;

% Run simulation and gather data
[utsol,usol,x,k]=Examle_1(dt,T,L,n,PDE_Plot);
dx=x(2)-x(1);

%% Now prepare the library data

% Prepare ux,uxx
ux=[];uxx=[];
ux_dum=[];uxx_dum=[];

% Take derivative using finite difference
for kk=1:size(usol,2)
    ux(kk,:)=CalDerivative(real(usol(:,kk)),dx,1);
    uxx(kk,:)=CalDerivative(real(usol(:,kk)),dx,2);
end

% Take derivative using FFT and compare the result of finite differene with
% this one.
ux_dum(1,:)=real(ifft((i*k).^1.*utsol(1,:)'));
uxx_dum(1,:)=real(ifft((i*k).^2.*utsol(1,:)'));

% Choose whether you want to plot
if PDE_Plot==1
    figure(1)
    clf
    hold on
    plot(real(usol(:,1)'),'linewidth',2.5,'color','blue')
    plot(ux(1,:),'linestyle','--','linewidth',2.5,'color','red')
    plot(uxx(1,:),'linestyle','-.','linewidth',2.5,'color','red')
    plot(ux_dum(1,:),'linestyle','--','linewidth',2.5,'color','green')
    plot(uxx_dum(1,:),'linestyle','-.','linewidth',2.5,'color','green')
    drawnow
end

% Reshape into a vector
Ux=reshape(ux',[],1);
Uxx=reshape(uxx',[],1);


%% Prepare the ut
ut=[];
for xpos=1:size(usol,1)
    dummy=real(usol(xpos,:));
    ut(:,xpos)=CalDerivative(dummy,dt,1);
    if PDE_Plot==1
        figure(1)
        clf
        hold on
        plot(ut(:,xpos))
    end
end

% Reshape
Ut=reshape(ut',[],1);

%% Prepare the u
u=real(usol)';

% Reshape
U=reshape(u',[],1);

%% Plot the result
if PDE_Plot==1
    figure(1)
    waterfall(x,tspan,real(usol'))
    colormap([0 0 0]);
    view(15,65)
    
    figure(2)
    waterfall(x,tspan,real(ux))
    colormap([0 1 0])
    view(15,65)
    
    figure(3)
    waterfall(x,tspan,real(uxx))
    colormap([0 0 1])
    view(15,65)
end

%% Now build your LHS and RHS library
syms U_sym Ut_sym Ux_sym Uxx_sym

LHS=[Ut Ut.*U Ut.*Ux Ut.*Ux.^2 Ut.*Ux.^3];
LHS_sym=[Ut_sym Ut_sym*U_sym Ut_sym*Ux_sym Ut_sym*Ux_sym^2 Ut_sym*Ux_sym^3];

RHS=[ones(size(U)) U Ux Uxx  Ux.^2 Uxx.^2  U.*Ux  U.*Uxx  U.^2.*Ux  U.^2.*Uxx ];
Symbol=[1,U_sym,Ux_sym,Uxx_sym,Ux_sym^2, Uxx_sym^2, U_sym*Ux_sym,  U_sym*Uxx_sym,  U_sym^2*Ux_sym,  U_sym^2*Uxx_sym];

RHS_sym=cell(1,10);
for ii=1:10
    RHS_sym{1,ii}=Symbol(1,ii);
end

Data_Length=size(U,1);
Random_Sequence=randperm(Data_Length)';

LHS=LHS(Random_Sequence,:);
RHS=RHS(Random_Sequence,:);

%% Now run SINDy-PI
% Define some parameter for sparese regression
lambda=0.3;
N=15;
disp=1;
NormalizeLib=1;

[Xi,ODEs] = sparsifyDynamics(RHS,LHS(:,2),LHS_sym(2),lambda,N,RHS_sym,disp,NormalizeLib);







