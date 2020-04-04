%% This function will calculate the prediction error between the discovered system and actual test data
% There are two methods to calculate the prediction error.
%
% Method one: This method suits for the one dimensional system only. The idea is
% using ODE45 to drive each point k-steps forward in time and calculate the
% prediction error between the system estimation and actual states. Here,
% the system's estimate is the value of derivative at k-steps.
%
% Method two: This method suits for one dimensional system and
% multi-dimensional system. We use immediate prediction of the system's
% derivative and calculate the prediction error of the system model's
% derivative output and actual derivative.
%
% Last Updated: 2019/06/04
% Coded By: K
%%
function [Score]=Get_Score(dData_test,Data_test,u,Control,tspan,Shuffle,name,Prediction_Steps,dt,size1,size2,method)
% Get the ODE simulation result
Noise=0;

% Start calculate
if method ==1
    % If the method one is selected, then we use the multi step prediction
    
    % Generate the function handel
    if u==0
        ODE_func=str2func(strcat('@(t,z)',name,'(t,z)'));
    end
    
    % Pass the test data as dummy variable
    Dummy=Data_test;
    
    if Prediction_Steps==0
        % This will result immediate prediction
        dData_Es=ODE_func(0,Dummy);
        Error=norm(dData_test-dData_Es);
        Dem=norm(dData_test);
        Score=Error/Dem;
    else
        dData_Es=ODE_func(0,Dummy);
        % Get the function estimation using RK45
        for pinpin=1:Prediction_Steps
            RK_k1=dt*dData_Es;
            RK_k2=dt*ODE_func(0,Dummy+0.5*RK_k1);
            RK_k3=dt*ODE_func(0,Dummy+0.5*RK_k2);
            RK_k4=dt*ODE_func(0,Dummy+RK_k3);
            Dummy=Dummy+(1/6)*(RK_k1+2*RK_k2+2*RK_k3+RK_k4);
            dData_Es=ODE_func(0,Dummy);
        end
        
        % Reshape
        Dum1=reshape(dData_Es,size1,size2);
        Dum2=reshape(dData_test,size1,size2);
        
        % Arrange all the data to the same time 
        Dum3=Dum1(1:end-Prediction_Steps,:);
        Dum4=Dum2(1+Prediction_Steps:end,:);
        
        % Calculate the prediction error
        Error=norm(reshape((Dum4-Dum3),[],1));
        Dem=norm(reshape(Dum4,[],1));
        
        % Calculate the norm
        Score=Error/Dem;
    end
else
    % Generate the function handel
    if u==0
        ODE_func=str2func(strcat('@(t,z)',name,'(t,z)'));
    end
   
    % Get the function estimation
    dData_Es=ODE_func(0,Data_test);
    Error=norm(dData_test-dData_Es);
    Dem=norm(dData_test);
    Score=Error/Dem;
end



