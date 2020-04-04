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
DictName_DL='Results/Result_DL_SINDy_Data_Length_Compare_State_6_LHS_Guess_2_V2';
mat_DL = dir(strcat(DictName_DL,'/*.mat'));

% Which LHS you are using
LHS_Guess=2;

% Determine the corretc vector that corresponding to the LHS you are
% choosing
if LHS_Guess==1
    Correct_Vector=LHS_Guess_1_Correct_Vector;
elseif LHS_Guess==2
    Correct_Vector=LHS_Guess_2_Correct_Vector;
end

%% Extract the result
% We read first file and extract some parameter information
for total_DL =1:1
    % Load the result from the folder
    load(strcat(DictName_DL,'/',mat_DL(total_DL).name));
    
    % Get the x axis, determine which percentage you used
    x_axis=percent_start:d_percent:percent_end;
    
    % Get the Lambada you have
    Lambda=(lam_start:lam_end)*d_lambda;
end


% Create a vector to store whether you got the structure correct
is_Right_DL=zeros(size(Lambda,2),size(x_axis,2));

% We read each result sperately
for total_DL = 1:length(mat_DL)
    % Load the result from the folder
    load(strcat(DictName_DL,'/',mat_DL(total_DL).name));
    
    % Loop through the percentage
    for i=1:size(Xi,1)
        % Loop through the lambda
        for j=1:size(Xi,2)
            % We now analyze the parameter vector Xi.
            Xi_dum=Xi{i,j};
            Xi_Active_dum=Xi_dum~=0;
            
            if Xi_Active_dum==Correct_Vector
                % Add one if the structure is correct
                is_Right_DL(j,i)=is_Right_DL(j,i)+1;
            end
        end
    end
end


%%
figure(1)
hold on
yValue=max(is_Right_DL)/total_DL*100;
xValue=x_axis*100;
XlabelTick=cell(size(xValue));
YlabelTick=cell(size(yValue));

plot(xValue,yValue,'linewidth',2,'color','blue')
scatter(xValue,yValue,50,'filled','MarkerFaceColor','r')

set(gca,'FontSize',30);
ax_x=cellstr(num2str(get(gca,'xtick')'))
ax_y=cellstr(num2str(get(gca,'ytick')'))
for i=1:size(ax_x,1)
    XlabelTick{1,i}=strcat(ax_x{i,1},'%');
end
for i=1:size(ax_y,1)
    YlabelTick{1,i}=strcat(ax_y{i,1},'%');
end

set(gca,'XTickLabel',XlabelTick,'YTickLabel',YlabelTick);
set(gcf,'Position',[100 100 600 400]);
set(gcf,'PaperPositionMode','auto');
grid on



