%% This function will generate the simulation data of an ODE function.
% Last Update: 2019/04/21
% Coded By: K

function [d_Data,Data]=Get_Sim_Data(ODE,state0,u,tspan,Noise,Control,Shuffle)
%% Get the size of the state and control
[N1,M1]=size(state0);
[N2,M2]=size(u);

%% Get simulation data by simulating the system using ODE113

% Determine the left hand side derivative
if Control==1
    y_list(1,:)=state0;
    d_y_list(1,:)=ODE(0,y_list(1,:),u(1,:));
    for i=2:length(u)
        [t_1,y_1] = ode113(@(t_1,y_1)ODE(t_1,y_1,u(i-1,:)),tspan(1,i-1:i),state0);
        y_list(i,:)=y_1(end,:);
        d_y_list(i,:)=ODE(0,y_list(i,:),u(i,:));
        state0=y_list(i,:)';
    end
else
    %Simulate the system ODE
    [t,y]=ode45(@(t,y)ODE(t,y),tspan,state0);
    y_list=y;
    % Get the derivative data
    d_y_list=ODE(0,y_list')';
end

%% Add some noise to the system
for i=1:N1
    Data(:,i)=y_list(:,i)+Noise*randn(size(y_list(:,i)));
end
%
for i=1:N1
    d_Data(:,i)=d_y_list(:,i)+Noise*randn(size(d_y_list(:,i)));
end

%% Shuffle the data
if Shuffle==1
    Sequence=randperm(size(Data,1));
    Data=Data(Sequence,:);
    d_Data=d_Data(Sequence,:);
end
