%% This function will plot the result of SINDy-PI peformance under different data length.
% It is meant to be used with YeastGlycolysis_DL_Auto_Test_Main_Part1.
% Coded By: K
% Last Updated: 2019/06/09
%% Clear all the variables
clc;clear all; close all;

%% Load the correct vector that SINDy-PI should discover
addpath('Datas')
addpath('Results')
load('Correct_Active_Term.mat')

%% Read the result from the specified folder

% Which LHS you are using
LHS_Guess=2;

if LHS_Guess==1
    DictName_DL='Results/Result_DL_SINDy_Data_Length_Compare_State_6_LHS_Guess_1_SwipeLength_Lambda_0.1';
    mat_DL = dir(strcat(DictName_DL,'/*.mat'));
elseif LHS_Guess==2
    DictName_DL='Results/Result_DL_SINDy_Data_Length_Compare_State_6_LHS_Guess_2_SwipeLength_Lambda_0.2';
    mat_DL = dir(strcat(DictName_DL,'/*.mat'));
end

% Determine the corretc vector that corresponding to the LHS you are
% choosing
if LHS_Guess==1
   Correct_Vector=LHS_Guess_1_Correct_Vector;
elseif LHS_Guess==2
    Correct_Vector=LHS_Guess_2_Correct_Vector;
end

%% Extract the result
% We read first file and extract some parameter information
for total_DL =1:length(mat_DL)
    % Load the result from the folder
    load(strcat(DictName_DL,'/',mat_DL(total_DL).name));
    
    Per_Vec(total_DL,1)=percent;
end

% Get the x axis, determine which percentage you used
x_axis=unique(Per_Vec);

% Create a vector to store how many time each percent is ran
Num_Per=zeros(size(x_axis));

% Create a vector to store whether you got the structure correct
is_Right_DL=zeros(size(x_axis));

% We read each result sperately
for total_DL = 1:length(mat_DL)
    % Load the result from the folder
    load(strcat(DictName_DL,'/',mat_DL(total_DL).name));
    
    % We now analyze the parameter vector Xi.
    Xi_dum=Xi;
    Xi_Active_dum=Xi_dum~=0;
    
    % Determine the index
    Index_Vec=x_axis==percent;
    [maxVal,maxPos]=max(Index_Vec);
    
    % Count how many times this percent shows up
    if Num_Per(maxPos,1)<=40
        Num_Per(maxPos,1)=Num_Per(maxPos,1)+1;
        
        if Xi_Active_dum==Correct_Vector
            % Add one if the structure is correct
            is_Right_DL(maxPos,1)=is_Right_DL(maxPos,1)+1;
        end
    end
end

%%
figure(1)
hold on
yValue=(is_Right_DL./Num_Per)*100;
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

% figure(2)
% plot(x_axis,Num_Per,'linewidth',2)
% ylim([0 40])