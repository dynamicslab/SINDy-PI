% Copyright 2015, All Rights Reserved
% Code by Steven L. Brunton
% Modified by Niall M. Mangan March 2016
% consider yeast glycolysis
% this file creates the data needed to recover each state variable
% once you have run this file, run the state-variable specific files:
% S1_yeast_glycolysis
% S2_yeast_glycolysis 
% S2_yeast_glycolysis_pareto (only this one does full pareto sweep and
%                               takes a very long time to run)
% S3_yeast_glycolysis (uses up to 3rd order polynomials in library)
% S4_yeast_glycolysis (uses up to 3rd order polynomials in library)
% S5_yeast_glycolysis (uses up to 3rd order polynomials in library)
% S6_yeast_glycolysis (Coming in later commit, requires more data, change sss value below)
% S7_yeast_glycolysis (uses up to 3rd order polynomials in library)
% Original Code: https://github.com/niallmm/iSINDy
%%
% This file is modified to generate the data.
% Modified by:K, 2019/07/16
%%
clear all; close all; clc; 
figpath = '../figures/';
addpath('./utils');
addpath('./bioutils');
changeplot;

plottag = 2;

%% generate Data

%size of system
n = 7;

% Integrate
dt = 0.1;
tspan = 0:dt:5; % time vector
options = odeset('RelTol',1e-7,'AbsTol',1e-7);

% % initial condition vector from paper
% A_3init = 2.475;
% N_2init = 0.077;
% S_1init = 1.187;
% S_2init = 0.193;
% S_3init = 0.050;
% S_4init = 0.115;
% S_5init = 0.077;
% Sinit = [S_1init S_2init S_3init S_4init S_5init A_3init N_2init];

rng(1); % seed the random number generator

%% This following settings is used to generate the baseline data.
% sss = 300; % 3*sss is the number of intial conditions
% Sinit = rand(sss,7);
% Sinit = [Sinit; 2*rand(sss,7)];
% Sinit = [Sinit; 3*rand(sss,7)];

%% This is used for the iSINDy (Needs more data)
% % need this many initial conditions to find the 6th state variable equation
sss = 450; %5*sss is the number of initial conditions.
Sinit = 0.01*rand(sss,7);
Sinit = 0.05*rand(sss,7);
Sinit = [Sinit; 0.1*rand(sss,7)];
Sinit = [Sinit; 1*rand(sss,7)];
Sinit = [Sinit; 2*rand(sss,7)];
Sinit = [Sinit; 3*rand(sss,7)];
Sinit = [Sinit; 4*rand(sss,7)];
Sinit = [Sinit; 5*rand(sss,7)];


measure = length(Sinit);

for ii = 1:measure
    %integrate for each initial conditions
    [t1,x1] = ode45(@yeastglycolysisNM, tspan, Sinit(ii,:));
    %store each instance
    tt(:,ii) = t1;
    x(:,:,ii) = x1;
end


%% add noise

eps = 0; % no error used here
xn = x + eps*randn(size(x));

% compute Derivative
% calculate exactly and add error
xt = [];dxt= []; t = [];

for ll =1:measure
    for ii=1:length(tspan)
        dxf(ii,:,ll) = yeastglycolysisNM(t,xn(ii,:, ll));
    end
    epsdt = 0;
    dxf = dxf+epsdt*randn(size(dxf));
    
    dxt = [dxt; dxf(:,:,ll)];
    xt = [xt; xn(:,:,ll)];
    t = [t; tt(:, ll)];
end

num2plot = 1:5*length(tspan);

if plottag >1
    % Plot data and derivatives
    figure(2)
    plot(t(num2plot,:) ,xt(num2plot,:),'o')
    hold on
    xlabel('time')
    ylabel('concentrations')
    drawnow
    
    figure(6)
    plot(t(num2plot,:), dxt(num2plot,:), '.')
    hold on
    xlabel('time')
    ylabel('derivative of concentrations w/ time')
    drawnow
end
