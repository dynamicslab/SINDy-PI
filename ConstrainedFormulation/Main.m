%% This file is the example of constrained formulation of SINDy-PI
% Coded By: K
% Last Update: 2019/07/09
%%
clc;close all;clear all;
%%
% Define the parameters for simulation
dt=0.01;T=3;

tspan=0:dt:T;

jx=0.6;Vmax=1.5;Km=0.3;

x0=1;

% Simulate
[t,x]=ode45(@(t,x)MMK_ODE(t,x,jx,Vmax,Km),tspan,x0);
dx=MMK_ODE(0,x',jx,Vmax,Km)';

% Now build library
Theta=[ones(size(x)) x x.^2 dx dx.*x dx.*x.^2];

[n,m]=size(Theta);
Diag_m=eye(m,m);
% Create the initial guess of the
C0=Diag_m(:)';

% Set the threshold
lambda=0.2;

% Determine how many iteration you need
N=5;

% Begin the optimization problem
C=C0;
tic
for iter=1:N
    % Run the optimization
    cvx_begin quiet
        variable xi(m,m)
            minimize(norm(Theta-Theta*xi));
        subject to
            diag(xi)==zeros(m,1);
            if iter>1
                xi(smallinds)==0;
            end
    cvx_end
    % Use thresholding
    fprintf('\n Iteration %d\n',iter)
    Xi=full(xi)
    smallinds = (abs(Xi)<lambda);
end
toc










