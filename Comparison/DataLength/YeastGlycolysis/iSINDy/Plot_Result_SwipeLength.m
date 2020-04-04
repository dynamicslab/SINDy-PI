%% This function will plot the result of implicit-SINDy peformance under different data length.
% Coded By: K
% Last Updated: 2019/06/09
%% Clear all the variables
clc;clear all; close all;

%% Load the correct vector that implicit-SINDy should discover
addpath('Results')
load('Correct_Xi.mat')
Correct_Pin=Xi_True~=0;

%% Read the result from the specified folder
DictName_im='Results/';
mat_DL = dir(strcat(DictName_im,'/*.mat'));

%% Extract the result
% We read first file and extract some parameter information
for total_Im =1:length(mat_DL)
    % Load the result from the folder
    load(strcat(DictName_im,'/',mat_DL(total_Im).name));
    
    Per_Vec(total_Im,1)=percent;
end

% Get the x axis, determine which percentage you used
x_axis=uniquetol(Per_Vec,0.001);

% Create a vector to store how many time each percent is ran
Num_Per=zeros(size(x_axis));

% Create a vector to store whether you got the structure correct
is_Right_Im=zeros(size(x_axis));

% We read each result sperately
for total_Im = 1:length(mat_DL)
    % Load the result from the folder
    load(strcat(DictName_im,'/',mat_DL(total_Im).name));
    mat_DL(total_Im).name
    % We now analyze the parameter vector Xi.
    Xi_dum=Xi1;
    Xi_Active_dum=Xi_dum~=0;
    
    % Determine the index
    Minus_Result=abs(x_axis-percent);
    for i=1:length(Minus_Result)
        if Minus_Result(i,1)<0.001
            Minus_Result(i,1)=1;
        else
            Minus_Result(i,1)=0;
        end
    end
    Index_Vec=Minus_Result;
    [maxVal,maxPos]=max(Index_Vec);
    
    Num_Per(maxPos,1)=Num_Per(maxPos,1)+1;
    
    if sum(abs(Xi_Active_dum-Correct_Pin))==0
        % Add one if the structure is correct
        is_Right_Im(maxPos,1)=is_Right_Im(maxPos,1)+1;
    end
    
end

%%
figure(1)
hold on
yValue=(is_Right_Im./Num_Per)*100;
xValue=x_axis*100;

xValue=x_axis*100;

d_step=0.01;
p=pchip(xValue,yValue,min(xValue):d_step:max(xValue));

plot(min(xValue):d_step:max(xValue),p,'linewidth',3,'color','black','linestyle','-')
scatter(xValue,yValue,200,'filled',...
    'MarkerFaceColor',0.4*[1 1 1],...
    'MarkerEdgeColor',[0 0 0])
for i=1:length(yValue)
    if yValue(i)==100
        break
    end
end

plot([xValue(i) xValue(i)],[0 100],'linewidth',3,'linestyle','-.','color','k')

% plot(xValue,yValue,'linewidth',2,'color','blue')
% scatter(xValue,yValue,50,'filled','MarkerFaceColor','r')
grid on
% title('Success Rate vs Data Usage','FontSize',18)
% xlabel('Data Usage','FontSize', 18)
% ylabel('Success Rate','FontSize', 18)
set(gcf,'Position',[100 100 600 400]);
set(gcf,'PaperPositionMode','auto');
set(gca,'FontSize',28);

ax_x=cellstr(num2str(get(gca,'xtick')'))
ax_y=cellstr(num2str(get(gca,'ytick')'))
for i=1:size(ax_x,1)
    XlabelTick{1,i}=strcat(ax_x{i,1},'%');
end
for i=1:size(ax_y,1)
    YlabelTick{1,i}=strcat(ax_y{i,1},'%');
end
set(gca,'XTickLabel',XlabelTick,'YTickLabel',YlabelTick);
box('on')



