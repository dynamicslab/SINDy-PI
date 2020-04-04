%% This file will generate the necessary simulation data needed for the comparison of
% noise sensitivity of iSINDy and SINDy-PI. The output of the file is
% training data. The training data contains noise. We will use three
% different methods to generate the training data.
%
% Method one: This method will fake the derivative data. Instead of
% approximating the derivative data using the noisy actual data, we will
% directly add Guassian noise to the perfect derivative data.
%
% Method two: This method will directly take the finite difference on the
% noise data and approximate the derivative.
%
% Method three: This method will take the derivative based on TVRegDiff
% function.
%
% Date: 2019/06/04
% Coded By: K

function [xt,dxt,iter_length]=Generate_Noisy_Data(noise_level,init_num,dt,T,Denoise_Method,Gap)
%% Simulate the budworm population growth and gather the simulation data

% Define the system parameters
jx=0.6;Vmax=1.5;Km=0.3;

% Determine the simulation time step and time span
tspan=0:dt:T;

% Define noise level and add gaussian noise to the data
noise=0;

% Define whehter you have control, if you have it, please define it
Control=0;u=0;

%Define whether you want to shuffel the final data
Shuffle=0;

%% Generate the training data
if Denoise_Method==1
    %Method one: We fake the derivative data by directly adding guassian
    %noise to the perfect derivative data.
    
    % Define a dummy parameter to store the calculation
    dxt=[];xt=[];
    % Define the length of the data. Note that we will only use the middel
    % 40% of the data. Thus the first and last 30% of the data will be
    % discarded.
    len=length(tspan);
    start=round(0.3*len);
    iter_length=length(start:len-start);
    
    % Now adding the noise to the data
    if noise_level==0
        % If the noise level is zero, we will directly use the perfect
        % data.
        parfor iter=1:init_num
            % Randomlize the initial condition
            x0(iter,1)=abs(((round(iter/Gap))+0.5)*rand);
            % Run the ODE file
            [dxt_dum,xt_dum]=Get_Sim_Data(@(t,y)MMK_ODE(t,y,jx,Vmax,Km),x0(iter,1),u,tspan,noise,Control,Shuffle);
            % Now stack the value of the training data
            dxt=[dxt;dxt_dum(start:end-start)];
            xt=[xt;xt_dum(start:end-start)];
        end
    else
        parfor iter=1:init_num
            % Randomlize the initial condition
            x0(iter,1)=abs(((round(iter/Gap))+0.5)*rand);
            % Run the ODE file
            [dxt_dum,xt_dum]=Get_Sim_Data(@(t,y)MMK_ODE(t,y,jx,Vmax,Km),x0(iter,1),u,tspan,noise,Control,Shuffle);
            % Add Guassian noise to perfect data
            xt_dum2=xt_dum+noise_level*randn(size(xt_dum));
            dxt_dum2=dxt_dum+noise_level*randn(size(xt_dum));
            % Now stack the value of the training data
            dxt=[dxt;dxt_dum2(start:len-start)];
            xt=[xt;xt_dum2(start:len-start)];
        end
    end

elseif Denoise_Method==2
    % Use method two: We will directly take derivative with respect to
    % nosie data and approximate the derivative.
    
    % Define a dummy parameter to store the calculation
    dxt=[];xt=[];
    % Define the length of the data. Note that we will only use the middel
    % 40% of the data. Thus the first and last 30% of the data will be
    % discarded.
    len=length(tspan);
    start=round(0.3*len);
    iter_length=length(start:len-start);
    
    if noise_level==0
        % If the noise level is zero, we will directly use the perfect
        % data.
        parfor iter=1:init_num
            % Randomlize the initial condition
            x0(iter,1)=abs(((round(iter/Gap))+0.5)*rand);
            % Run the ODE file
            [dxt_dum,xt_dum]=Get_Sim_Data(@(t,y)MMK_ODE(t,y,jx,Vmax,Km),x0(iter,1),u,tspan,noise,Control,Shuffle);
            % Now stack the value of the training data
            dxt=[dxt;dxt_dum(start:end-start)];
            xt=[xt;xt_dum(start:end-start)];
        end
    else
        parfor iter=1:init_num
            % Randomlize the initial condition
            x0(iter,1)=abs(((round(iter/Gap))+0.5)*rand);
            % Run the ODE file
            [dxt_dum,xt_dum]=Get_Sim_Data(@(t,y)MMK_ODE(t,y,jx,Vmax,Km),x0(iter,1),u,tspan,noise,Control,Shuffle);
            % Add Guassian noise to the x
            xt_dum2=xt_dum+noise_level*randn(size(xt_dum));
            % Take finite difference on the x_noise. This could be achieved
            % by using the gradient command.
            dxt_dum=gradient(xt_dum2,dt);
            % Now stack the value of the training data
            dxt=[dxt;dxt_dum(start:len-start)];
            xt=[xt;xt_dum2(start:len-start)];
        end
    end
    
elseif Denoise_Method==3
    % Use method three: We will use the TVRegDiff method to calculate the
    % deravtive.
    
    % Define the parameter used for the TVRegDiff. The general idea is
    % using higher smoothing factor as the noise level increases.
    if noise_level <= 1e-3
        Alpha=1;
    elseif noise_level <= 0.002
        Alpha=5e1;
    elseif noise_level<=0.004
        Alpha=1e3;
    elseif noise_level<0.008
        Alpha=1.5e3;
    elseif noise_level<0.01
        Alpha=2e3;
    elseif noise_level<=0.01
        Alpha=3e3;
    elseif noise_level<=0.02
        Alpha=3.1e3;
    elseif noise_level<=0.03
        Alpha=3.5e3;
    elseif noise_level<=0.04
        Alpha=4e3;
    elseif noise_level<=0.05
        Alpha=5.2e3;
    elseif noise_level<=0.06
        Alpha=5.5e3;
    elseif noise_level<=0.07
        Alpha=8e3;
    elseif noise_level<=0.08
        Alpha=9e4;
    elseif noise_level<=0.09
        Alpha=9.5e4;
    elseif noise_level<=0.1
        Alpha=1e4;
    elseif noise_level<=0.13
        Alpha=1.2e4;
    elseif noise_level<=0.16
        Alpha=1.6e4;
    elseif noise_level<=0.19
        Alpha=2e4;
    elseif noise_level<=0.6
        Alpha=5e4;
    else
        Alpha=1e6;
    end
    
    % Define parameters for the TVRegDiff function
    plotflg=0;diagflag=0;
    
    % Define a dummy parameter to store the calculation
    dxt=[];xt=[];
    % Define the length of the data. Note that we will only use the middel
    % 40% of the data. Thus the first and last 30% of the data will be
    % discarded.
    len=length(tspan);
    start=round(0.3*len);
    iter_length=length(start:len-start);
    
    if noise_level==0
        % If the noise level is zero, we will directly use the perfect
        % data.
        parfor iter=1:init_num
            % Randomlize the initial condition
            x0(iter,1)=abs(((round(iter/Gap))+0.5)*rand);
            % Run the ODE file
            [dxt_dum,xt_dum]=Get_Sim_Data(@(t,y)MMK_ODE(t,y,jx,Vmax,Km),x0(iter,1),u,tspan,noise,Control,Shuffle);
            % Now stack the value of the training data
            dxt=[dxt;dxt_dum(start:end-start)];
            xt=[xt;xt_dum(start:end-start)];
        end
    else
        parfor iter=1:init_num
            % Randomlize the initial condition
            x0(iter,1)=abs(((round(iter/Gap))+0.5)*rand);
            % Run the ODE file
            [dxt_dum,xt_dum]=Get_Sim_Data(@(t,y)MMK_ODE(t,y,jx,Vmax,Km),x0(iter,1),u,tspan,noise,Control,Shuffle);
            % Add Guassian noise to the x
            xt_dum2=xt_dum+noise_level*randn(size(xt_dum));
            % Calculate the derivative using TVRegDiff
            dxt_dum2=TVRegDiff( xt_dum2, 10, Alpha, [], 'small', 1e8, dt, plotflg, diagflag );
            % Now stack the value of the training data
            dxt=[dxt;dxt_dum2(start:len-start)];
            xt=[xt;xt_dum2(start:len-start)];
        end
    end
end







