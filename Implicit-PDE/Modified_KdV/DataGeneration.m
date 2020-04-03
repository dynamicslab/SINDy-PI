%% This file is used to generate the simulation data of the Schwarzian-Kroteweg-de Vires equation.
% Coded By: K
% Last Updated: 2019/06/17
%% Clear all
clc;clear;close all;

%% Create a folder
addpath('Functions')
%
[fld_status, fld_msg, fld_msgID]=mkdir('Datas');
addpath('Datas')
%
[fld_status, fld_msg, fld_msgID]=mkdir('TempFunctions');
addpath('TempFunctions')
%% Define the parameters for the SKdV equation
% Define the time span
dt=0.01;T=20;
tspan=0:dt:T;

% Define the spatial domain
L=50; 

% Define the discretization point
n=256; 

% Get the index of each point
x2=linspace(-L/2,L/2,n+1); 
x=x2(1:n); 
dx=x(2)-x(1);

% Get the frequency
k=(2*pi/L)*[0:n/2-1 -n/2:-1].';
k2=fftshift(k);

% Get the initial value
a1=0.5; 
u1=(a1/2*(sech(0.5*sqrt(1)*(x+10))).^2 ).'; 

a2=0.5; 
u2=(a2/2*(sech(0.5*sqrt(1)*(x-10))).^2 ).'; 

u0=u1+u2;


% Transfer the initial point to spectram domain
ut0=fft(u0); 

% Define parameter for system
% Omega: This term should be positive, and it will determine the
% dissapation of the system
Omega = 0.1; 

alpha = 1;

% g0: This term should be positive, and it determines the gain of the
% system.
g0 =[0;0.0001;0.001;0.01;0.1;1]; 

%% Loop through the g0
for kkk=1:size(g0,1)
    % Solve for the PDE
    tic
    opts = odeset('RelTol',1e-12,'AbsTol',1e-13);
    [time,utsol]=ode45(@(time,u0)mkdv_rhs(time,u0,k,x,Omega,g0(kkk,1),alpha),tspan,ut0,opts);
    toc
    
    % Back to time domain
    for j=1:length(tspan)
        usol(:,j)=ifft(utsol(j,1:n).');
    end
    
    u=zeros(size(usol));
    ut=zeros(size(usol));
    ux=zeros(size(usol));
    uxx=zeros(size(usol));
    uxxx=zeros(size(usol));
    
    % Use the previous utsol to get the derivative
    for j=1:length(tspan)
        [~,ut(:,j),u(:,j),ux(:,j),uxx(:,j),uxxx(:,j)]=mkdv_rhs(0,utsol(j,1:n)',k,x,Omega,g0(kkk,1),alpha);
    end
    
    u=u';ut=ut';ux=ux';uxx=uxx';uxxx=uxxx';
    
    name=strcat('Datas/','SimulationData_V_',num2str(kkk),'.mat');
    save(name,'u','ut','ux','uxx','uxxx','x','tspan','dt','T','L','n');
        
end

%% Plot the result
figure(1)
waterfall(x,tspan,real(u))
colormap([0 0 0]);
view(42,55)

figure(2)
surfl(x,tspan,real(u))
colormap(gray)
shading interp
view(42,55)

figure(3)
pcolor(x,tspan,real(u))
shading interp




