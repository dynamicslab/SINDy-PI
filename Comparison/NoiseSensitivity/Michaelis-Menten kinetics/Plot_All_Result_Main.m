%% This function will plot whether SINDy-PI and implicit-SINDy finds out the correct structure.
% Coded By: K
% Last Updated: 2019/09/11
% Modified: The way to calculate the parameter error.
%% Clear all the variables
clc;clear all; close all;
%%
% Choose which denoise method you want:
% 1 is adding Guassian noise to actual derivative, 2 is direct
% finite difference, 3 is TVRegDiff.
addpath('Functions')
DictName_DL='Result_DL_SINDy_Derivative_Method2_PridictionStep_0_IniNum_2400';
DictName_iS='Result_i_SINDy_Derivative_Method2_PridictionStep_0_IniNum_2400';
mat_DL = dir(strcat(DictName_DL,'/*.mat'));
mat_iS = dir(strcat(DictName_iS,'/*.mat'));

%%
% Prepare the noise level vector for plot
percent_start=1;
percent_end=24;
percent=0;
for percent_iter=percent_start:percent_end
    % Set the pin one step forward
    percent=percent+1;
    
    % Set the noise level
    noise_level(percent,1)=Determine_Noise_Level(percent_iter);
end

%%
% SINDy-PI Correct Result Rate: This rate is the rate that SINDy-PI Generates the correct structure.
% It does not mean that the correct structure is been selected.

% We first determine whether the SINDy-PI yields the correct structural
% result.
Correct_Structure_Count_DL=zeros(size(noise_level,1),length(mat_DL));

% Define the corresponding sparse vector. For the SINDy-PI, when the LHS
% choice is dx, the correct sparse vector we get should be Coffs_LHS1. When
% the LHS choice is dx*x, we could have two possible sparse vector,
% Coffs_LHS2 and Coffs_LHS2_2. WHen the LHS choice is dx*sin(x), it is not
% possible to get the correct structure due the sin(x). The choice 3 is not
% a term in the actual equation of the MMK.
Coffs_LHS1=[0.6;-3;0;0;0;-10/3;0;0;0];
Coffs_LHS2=[0.18;-0.9;0;0;0;-0.3;0;0;0];
Coffs_LHS2_2=[0;0.6;-3;0;0;0;-10/3;0;0];

% We now exame all the SINDy-PI generated result, and see whether it
% generates the correct structure.
for total_DL = 1:length(mat_DL)
    % Load the result from the folder
    load(strcat(DictName_DL,'/',mat_DL(total_DL).name));
    
    for pin1=1:size(Xi,1)
        % The first loop, pin1, will loop through different noise level.
        for pin2=1:size(Xi,3)-1
            % The second loop, pin2, will loop through different LHS
            % choice. We only care about LHS choice one and two, thus we use size(Xi,3)-1.
            for pin3=1:size(Xi,4)
                % The third loop, pin3, will loop through the different
                % sparsity parameter.
                
                % Extract the identified Xi
                Identified_Xi=Xi{pin1,1,pin2,pin3};
                Identified_Ind=Identified_Xi~=0;
                
                % Now determine whether it is correct or not.
                if pin2==1
                    % This means the first LHS choice
                    Ground_Truth=Coffs_LHS1~=0;
                    
                    % If the identified index matches the actual true
                    % index, we set counter to 1.
                    if abs(sum(Identified_Ind-Ground_Truth))==0
                        Correct_Structure_Count_DL(pin1,total_DL)=1;
                    end
                elseif pin2==2
                    % This means the second LHS choice
                    Ground_Truth1=Coffs_LHS2~=0;
                    Ground_Truth2=Coffs_LHS2_2~=0;
                    
                    % If the identified index matches the actual true
                    % index, we set counter to 1.
                    if abs(sum(Identified_Ind-Ground_Truth1))==0 || abs(sum(Identified_Ind-Ground_Truth2))==0
                        Correct_Structure_Count_DL(pin1,total_DL)=1;
                    end
                end
            end
        end
    end
end

Correct_Result_Rate_DL=(sum(Correct_Structure_Count_DL')/total_DL)*100;

%%
% Implicit-SINDy Correct Result Rate: This rate is the rate that Implicit-SINDy Generates the correct structure.
% It does not mean the rate that correct structure is been selected.

% We first determine whether the Implicit-SINDy yields the correct structural
% result.
Correct_Structure_Count_Im=zeros(size(noise_level,1),length(mat_iS));

% Define the corresponding sparse vector. For the Im-SINDy, there are four
% possible correct solution.
Coffs_Ind1=[1;1;0;0;0;1;1;0;0;0];
Coffs_Ind2=[0;1;1;0;0;0;1;1;0;0];
Coffs_Ind3=[0;0;1;1;0;0;0;1;1;0];
Coffs_Ind4=[0;0;0;1;1;0;0;0;1;1];

% We now exame all the SINDy-PI generated result, and see whether it
% generates the correct structure.
for total_iS = 1:length(mat_iS)
    % Load the result from the folder
    load(strcat(DictName_iS,'/',mat_iS(total_iS).name));
    
    for pin1=1:size(Xi_ns,1)
        % The first loop, pin1, will loop through different noise level
        Identified_Im=Xi_ns{pin1,1};
        
        for pin2=1:size(Identified_Im,2)
            % pin2 will loop through all the possible sparse vector
            % identified by implicit-SINDy
            Possible_Im=Identified_Im(:,pin2);
            Identified_Ind_Im=Possible_Im~=0;
            
            if sum(abs(Identified_Ind_Im-Coffs_Ind1))==0
                Correct_Structure_Count_Im(pin1,total_iS)=1;
            elseif sum(abs(Identified_Ind_Im-Coffs_Ind2))==0
                Correct_Structure_Count_Im(pin1,total_iS)=1;
            elseif sum(abs(Identified_Ind_Im-Coffs_Ind3))==0
                Correct_Structure_Count_Im(pin1,total_iS)=1;
            elseif sum(abs(Identified_Ind_Im-Coffs_Ind4))==0
                Correct_Structure_Count_Im(pin1,total_iS)=1;
            end
        end
    end
end

Correct_Result_Rate_Im=(sum(Correct_Structure_Count_Im')/total_iS)*100;


%%
% SINDy-PI correct identification rate of the selected model(the rate that
% the model with lowest prediction error has the correct structure).

% We first extract the result of L2 prediction error
for total_DL = 1:length(mat_DL)
    % Load the result from the folder
    load(strcat(DictName_DL,'/',mat_DL(total_DL).name));
    
    % Store the minumum L2 prediction error of the SINDy-PI
    % Index1 will tell you which lambda generates the best model.
    % Index2 will tell you which left hand side generates the best model.
    % min_DLS_Val stores the minimum L2 prediction error of each run.
    for pin1=1:size(Score_DLS,1)
        for pin3=1:size(Score_DLS,3)
            [min_DLS_Val_LHS(pin1,pin3),Index1(total_DL,pin1,pin3)]=min(Score_DLS(pin1,1,pin3,:));
            [AA,Index]=min(Score_DLS(pin1,1,pin3,:))
        end
        [min_DLS_Val(total_DL,pin1),Index2(total_DL,pin1)]=min(min_DLS_Val_LHS(pin1,:));
    end
end

Total_Selection_LHS1=zeros(9,percent_end-percent_start+1);
Total_Selection_LHS2=zeros(9,percent_end-percent_start+1);
is_Right_DL=zeros(percent_end-percent_start+1,1);
LHS3_Count=0;
for total_DL = 1:length(mat_DL)
    % Load the result from the folder
    load(strcat(DictName_DL,'/',mat_DL(total_DL).name));
    
    for pin1=1:size(Score_DLS,1)
        Which_LHS=Index2(total_DL,pin1);
        Which_Lambda=Index1(total_DL,pin1,Which_LHS);
        Xi_Es_DL=Xi{pin1,1,Which_LHS,Which_Lambda};
        
        % For SINDy-PI, according to which LHS generates the best model, we
        % will have different corresponding correct parameters. Thus, the
        % calculation of differences between the parameter index norm will
        % differ.
        if Which_LHS==1
            Index_Norm_LHS1=Coffs_LHS1~=0;
            Index_Norm_Xi_Es_DL=Xi_Es_DL~=0;
            Index_Norm_Diff_DL(total_DL,pin1)=sum(abs(Index_Norm_LHS1-Index_Norm_Xi_Es_DL));
            Index_Norm_DL(total_DL,pin1)=sum(abs(Index_Norm_Xi_Es_DL));
            Coffs_Error_DL(total_DL,pin1)=norm(Coffs_LHS1-Xi_Es_DL,1);
            
            % If the structure is correct, we add one and count how many
            % times we get it right.
            if sum(abs(Index_Norm_Xi_Es_DL-Index_Norm_LHS1))==0
                is_Right_DL(pin1,1)=is_Right_DL(pin1,1)+1;
            else
                % When the discovered structure is not right, we calculate
                % which term is wrongly selected
                Total_Selection_LHS1(:,pin1)=Total_Selection_LHS1(:,pin1)+abs(Index_Norm_Xi_Es_DL-Index_Norm_LHS1);
            end
            
        elseif Which_LHS==2
            % If the left hand side is x*d(x), then there are two possible
            % solutions, we have to calculate them sperately.
            Index_Norm_Xi_Es_DL=Xi_Es_DL~=0;
            
            % If the structure is correct, we add one and count how many
            % times we get it right.
            if sum(abs(Index_Norm_Xi_Es_DL-(Coffs_LHS2~=0)))==0
                % Possible solution one:
                is_Right_DL(pin1,1)=is_Right_DL(pin1,1)+1;
                Index_Norm_DL(total_DL,pin1)=sum(abs(Index_Norm_Xi_Es_DL));
                Index_Norm_LHS2=Coffs_LHS2~=0;
                Index_Norm_Diff_DL(total_DL,pin1)=sum(abs(Index_Norm_LHS2-Index_Norm_Xi_Es_DL));
                Coffs_Error_DL(total_DL,pin1)=norm(Coffs_LHS2-Xi_Es_DL,1);
            elseif sum(abs(Index_Norm_Xi_Es_DL-(Coffs_LHS2_2~=0)))==0
                % Possible solution two:
                is_Right_DL(pin1,1)=is_Right_DL(pin1,1)+1;
                Index_Norm_DL(total_DL,pin1)=sum(abs(Index_Norm_Xi_Es_DL));
                Index_Norm_LHS2=Coffs_LHS2_2~=0;
                Index_Norm_Diff_DL(total_DL,pin1)=sum(abs(Index_Norm_LHS2-Index_Norm_Xi_Es_DL));
                Coffs_Error_DL(total_DL,pin1)=norm(Coffs_LHS2_2-Xi_Es_DL,1);
            else
                % The structure discovered is not correct
                Index_Norm_DL(total_DL,pin1)=sum(abs(Index_Norm_Xi_Es_DL));
                Index_Norm_LHS2=Coffs_LHS2~=0;
                Index_Norm_Diff_DL(total_DL,pin1)=sum(abs(Index_Norm_LHS2-Index_Norm_Xi_Es_DL));
                Coffs_Error_DL(total_DL,pin1)=norm(Coffs_LHS2-Xi_Es_DL,1);
                
                % When the discovered structure is not right, we calculate
                % which term is wrongly selected
                Total_Selection_LHS2(:,pin1)=Total_Selection_LHS2(:,pin1)+abs(Index_Norm_LHS2-Index_Norm_Xi_Es_DL);
                
            end
            
        elseif Which_LHS==3
            % Count how many times LHS 3 is selected
            LHS3_Count=LHS3_Count+1;
        end
    end
end



%%
% Im-SINDy correct identification rate of the selected model(the rate that
% the model with lowest prediction error has the correct structure).

% We first extract the result of L2 prediction error
is_Right_iS=zeros(percent_end-percent_start+1,1);

for total_iS = 1:length(mat_iS)
    load(strcat(DictName_iS,'/',mat_iS(total_iS).name));
    
    % Store the minimum of implicit SINDy prediction norm, also its
    % corresponding index.
    for pin1=1:size(Score_iS,1)
        [min_iS_Val(total_iS,pin1),Index1(total_iS,pin1)]=min(Score_iS(pin1,:));
    end
    
    % Test whether the correct structure is identified
    for pin1=1:size(Score_iS,1)
        % Extract the possible null space vector for each noise level
        Xi_dum=Xi_ns{pin1,1};
        % Extract the best parameter vector that produces the minimum prediction error
        Xi_Es_iS=Xi_dum(:,Index1(total_iS,pin1));
        % Get the activated terms in this vector
        Coffs_dum=Xi_Es_iS~=0;
        % The implicit-SINDy might generate several solutions that results
        % the correct structure after the symbolic simplification. To
        % calculate the parameter error we need to consider each case
        % sperately.
        
        % First calculate the index norm
        Index_Norm_iS(total_iS,pin1)=sum(abs(Coffs_dum));
        
        % Then we determine whether the correct structure is identified.
        % For the implicit-SINDy, there are many  different possible null
        % space vector that will generates the correct structure. We need
        % to test them one by one.
        
        % The first possible solution.
        if sum(abs(Coffs_dum-[1;1;0;0;0;1;1;0;0;0]))==0
            is_Right_iS(pin1,1)=is_Right_iS(pin1,1)+1;
            Coffs_norm=Xi_Es_iS/Xi_Es_iS(1);
            Coffs_Error_iS(total_iS,pin1)=norm(Coffs_norm-[-0.18;0.9;0;0;0;0.3;1;0;0;0]/(-0.18),1);
            Index_Norm_Diff_iS(total_iS,pin1)=0;
            % The second possible solution.
        elseif sum(abs(Coffs_dum-[0;1;1;0;0;0;1;1;0;0]))==0
            is_Right_iS(pin1,1)=is_Right_iS(pin1,1)+1;
            Coffs_norm=Xi_Es_iS/Xi_Es_iS(2);
            Coffs_Error_iS(total_iS,pin1)=norm(Coffs_norm-[0;-0.18;0.9;0;0;0;0.3;1;0;0]/(-0.18),1);
            Index_Norm_Diff_iS(total_iS,pin1)=0;
            % The third possible solution.
        elseif sum(abs(Coffs_dum-[0;0;1;1;0;0;0;1;1;0]))==0
            is_Right_iS(pin1,1)=is_Right_iS(pin1,1)+1;
            Coffs_norm=Xi_Es_iS/Xi_Es_iS(3);
            Coffs_Error_iS(total_iS,pin1)=norm(Coffs_norm-[0;0;-0.18;0.9;0;0;0;0.3;1;0]/(-0.18),1);
            Index_Norm_Diff_iS(total_iS,pin1)=0;
            % The fourth possible solution.
        elseif sum(abs(Coffs_dum-[0;0;0;1;1;0;0;0;1;1]))==0
            is_Right_iS(pin1,1)=is_Right_iS(pin1,1)+1;
            Coffs_norm=Xi_Es_iS/Xi_Es_iS(4);
            Coffs_Error_iS(total_iS,pin1)=norm(Coffs_norm-[0;0;0;-0.18;0.9;0;0;0;0.3;1]/(-0.18),1);
            Index_Norm_Diff_iS(total_iS,pin1)=0;
        else
            % Else, none of the solution is correct. We will calculate the
            % parameter error and index norm difference using the first
            % possible solutions.
            
            % We first try to figure out which solution the implicit-SINDy
            % is trying to approximate.
            Vec1(1,1)=Xi_Es_iS(1,1);
            Vec1(2,1)=Xi_Es_iS(2,1);
            Vec1(3,1)=Xi_Es_iS(6,1);
            Vec1(4,1)=Xi_Es_iS(7,1);
            %
            Vec2(1,1)=Xi_Es_iS(2,1);
            Vec2(2,1)=Xi_Es_iS(3,1);
            Vec2(3,1)=Xi_Es_iS(7,1);
            Vec2(4,1)=Xi_Es_iS(8,1);
            %
            Vec3(1,1)=Xi_Es_iS(3,1);
            Vec3(2,1)=Xi_Es_iS(4,1);
            Vec3(3,1)=Xi_Es_iS(8,1);
            Vec3(4,1)=Xi_Es_iS(9,1);
            %
            Vec4(1,1)=Xi_Es_iS(4,1);
            Vec4(2,1)=Xi_Es_iS(5,1);
            Vec4(3,1)=Xi_Es_iS(9,1);
            Vec4(4,1)=Xi_Es_iS(10,1);
            % Calculate the corresponding error
            Error_Dum(1,1)=norm(Vec1/Vec1(1,1)-[-0.18;0.9;0.3;1]/(-0.18),1)
            Error_Dum(2,1)=norm(Vec2/Vec2(1,1)-[-0.18;0.9;0.3;1]/(-0.18),1);
            Error_Dum(3,1)=norm(Vec3/Vec3(1,1)-[-0.18;0.9;0.3;1]/(-0.18),1);
            Error_Dum(4,1)=norm(Vec4/Vec4(1,1)-[-0.18;0.9;0.3;1]/(-0.18),1);
            % Get the minimal error and coresponding index
            [Min_Error_Dum,Min_Error_Dum_Index]=min(Error_Dum);
            % Now calculate the parameter error and index norm difference
            if Min_Error_Dum_Index==1
                Coffs_norm=Xi_Es_iS;
                Coffs_Error_iS(total_iS,pin1)=norm(Coffs_norm-[-0.18;0.9;0;0;0;0.3;1;0;0;0],1);
                Index_Norm_Diff_iS(total_iS,pin1)=sum(abs([1;1;0;0;0;1;1;0;0;0]-Coffs_dum));
                % The second possible solution.
            elseif Min_Error_Dum_Index==2
                Coffs_norm=Xi_Es_iS;
                Coffs_Error_iS(total_iS,pin1)=norm(Coffs_norm-[0;-0.18;0.9;0;0;0;0.3;1;0;0],1);
                Index_Norm_Diff_iS(total_iS,pin1)=sum(abs([0;1;1;0;0;0;1;1;0;0]-Coffs_dum));
                % The third possible solution.
            elseif Min_Error_Dum_Index==3
                Coffs_norm=Xi_Es_iS;
                Coffs_Error_iS(total_iS,pin1)=norm(Coffs_norm-[0;0;-0.18;0.9;0;0;0;0.3;1;0],1);
                Index_Norm_Diff_iS(total_iS,pin1)=sum(abs([0;0;1;1;0;0;0;1;1;0]-Coffs_dum));
                % The fourth possible solution.
            elseif Min_Error_Dum_Index==4
                Coffs_norm=Xi_Es_iS;
                Coffs_Error_iS(total_iS,pin1)=norm(Coffs_norm-[0;0;0;-0.18;0.9;0;0;0;0.3;1],1);
                Index_Norm_Diff_iS(total_iS,pin1)=sum(abs([0;0;0;1;1;0;0;0;1;1]-Coffs_dum));
            end
        end
        
    end
    
end

%% Plot the result for current iteration
close all
%% This figure will plot the success rate of the SINDy-PI and
% implicit-SINDy.
figure(1)
noise_level(1,1)=1e-8;
hold on
% This plots the rate that SINDy-PI best model does not have the correct structure
% pp1=plot(noise_level,100-(is_Right_DL/length(mat_DL))*100,'linewidth',2,'linestyle','-.','color','blue')
% scatter(noise_level,100-(is_Right_DL/length(mat_DL))*100,50,'blue')

% This plots the rate that SINDy-PI does not generate the correct structure
pp3=plot(noise_level,100-Correct_Result_Rate_DL,'linewidth',5,'linestyle','-','Color','blue')
scatter(noise_level,100-Correct_Result_Rate_DL,80,'filled','MarkerFaceColor','b')

% This plots the rate that im-SINDy best model does not have the correct structure
% pp2=plot(noise_level,100-is_Right_iS/total_iS*100,'linewidth',2,'linestyle','-.','color','red')
% scatter(noise_level,100-is_Right_iS/total_iS*100,50,'r')

% This plots the rate that im-SINDy does not generate the correct structure
pp3=plot(noise_level,100-Correct_Result_Rate_Im,'linewidth',5,'linestyle','-','Color','red')
scatter(noise_level,100-Correct_Result_Rate_Im,80,'filled','MarkerFaceColor','red')

xlim([1e-8 1])
drawnow
grid on
% title('Success Rate','FontSize',18)
% xlabel('Noise Level $\sigma$','FontSize', 18)
% ylabel('Success Rate','FontSize', 18)
% legend([pp1,pp2],'SINDy-PI','im-SINDy','location','NorthWest')
set(gca,'FontSize',28);
set(gcf,'Position',[100 100 2400 150]);
%set(gcf,'Position',[100 100 3400 250]);
set(gcf,'PaperPositionMode','auto');
set(gca, 'XScale', 'log')
box('on')
ax_y=cellstr(num2str(get(gca,'ytick')'))
%ax_x=cellstr(num2str(get(gca,'xtick')'))
for i=1:size(ax_y,1)
    YlabelTick{1,i}=strcat(ax_y{i,1},'%');
end
% Set the first point as zero
% for i=1:size(ax_x,1)
%     XlabelTick{1,i}=ax_x{i,1};
% end
XlabelTick=get(gca,'XTickLabel');
XlabelTick{1,1}=strcat('0');

set(gca,'YTickLabel',YlabelTick);
set(gca,'XTickLabel',XlabelTick);
%print('-depsc2', '-loose', 'Figures/DL-SINDy_vs_iSINDy_SuccessRate.eps');

%%
% This figure will plot the l2 prediction error of the SINDy-PI and
% implicit-SINDy.
figure(2)
hold on
plot1=plot(noise_level,mean(min_DLS_Val),'linewidth',5,'linestyle','-','Color','blue')
scatter(noise_level,mean(min_DLS_Val),80,'filled','MarkerFaceColor','b')
%
plot2=plot(noise_level,mean(min_iS_Val),'linewidth',5,'linestyle','-','Color','red')
scatter(noise_level,mean(min_iS_Val),80,'filled','MarkerFaceColor','r')
drawnow
grid on
% title('Peformance Comparison','FontSize',18)
% xlabel('Noise Level $\sigma$','FontSize', 18)
% ylabel('Prediction Error ($l_2$)','FontSize', 18)
% legend([plot1 plot2],'SINDy-PI','im-SINDy','location','NorthWest')
xlim([1e-8 1])

ax_y=cellstr(num2str(get(gca,'ytick')'))
ax_x=cellstr(num2str(get(gca,'xtick')'))

YlabelTick={1e-12,1e-6,1};

set(gca,'FontSize',28);
set(gcf,'Position',[100 100 2400 150]);
set(gcf,'PaperPositionMode','auto');
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
set(gca,'XTickLabel',[]);
set(gca,'YTick',[10^(-12) 10^(-6) 10^0]);
box('on')

%print('-depsc2', '-loose', 'Figures/DL-SINDy_vs_iSINDy_L2_Norm.eps');

%%
% This figure will plot the average l1 parameter error of the SINDy-PI and
% implicit-SINDy.
figure(3)
hold on
% Here we exclude the result which cause the parameter error becomes NaN
for ii=1:size(Coffs_Error_DL,2)
    pinpin=0;
    dummy_Sum=0;
    for jj=1:size(Coffs_Error_DL,1)
        if isnan(Coffs_Error_DL(jj,ii))==0
            pinpin=pinpin+1;
            dummy_Sum=dummy_Sum+Coffs_Error_DL(jj,ii);
        end
    end
    mean_Coffs_Error_DL(ii,1)=dummy_Sum/pinpin;
end

% Plot the parameter error of the SINDy-PI
ppp1=plot(noise_level,mean_Coffs_Error_DL,'linewidth',5,'linestyle','-','color','blue')
scatter(noise_level,mean_Coffs_Error_DL,80,'filled','MarkerFaceColor','b')
% Here we exclude the result which cause the parameter error becomes NaN
for ii=1:size(Coffs_Error_iS,2)
    pinpin=0;
    dummy_Sum=0;
    for jj=1:size(Coffs_Error_iS,1)
        if isnan(Coffs_Error_iS(jj,ii))==0
            pinpin=pinpin+1;
            dummy_Sum=dummy_Sum+Coffs_Error_iS(jj,ii);
        end
    end
    mean_Coffs_Error_iS(ii,1)=dummy_Sum/pinpin;
end


ppp2=plot(noise_level,mean_Coffs_Error_iS,'linewidth',5,'linestyle','-','color','red')
scatter(noise_level,mean_Coffs_Error_iS,80,'filled','MarkerFaceColor','r')
drawnow
grid on
% legend([ppp1 ppp2],'SINDy-PI','im-SINDy')
% title('Average parameter error','FontSize',18)
% xlabel('Noise Level $\sigma$','FontSize', 18)
% ylabel('Average Parameter Error ($l_1$)','FontSize', 18)
xlim([1e-8 1])
%ylim([1e-17 100])
set(gca,'FontSize',28);
set(gcf,'Position',[100 100 2400 150]);
set(gcf,'PaperPositionMode','auto');
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
set(gca,'XTickLabel',[]);
set(gca,'YTick',[10^(-12) 10^(-6) 10^0]);
box('on')

%%
% We will plot the index norm using boxplot
% Prepare the lable
noise_level(1,1)=0;
for qq=1:size(noise_level,1)
    Xlable{1,qq}=num2str(noise_level(qq,1))
end

figure(4)
hold on
ppp1=violin(Index_Norm_DL,'facecolor','blue','edgecolor','k',...
    'bw',0.4,'mc',[],'medc','k','facealpha','0.5')
plot(0:25,3*ones(1,size(noise_level,1)+2),'--','linewidth',2)
rectangle('Position',[0 2 25 2],'EdgeColor','none','FaceColor',[0 1 0 0.4])

chi=get(gca, 'Children')
set(gca, 'Children',flipud(chi))
legend('off')
grid on
set(gca,'FontSize',20);
set(gcf,'Position',[100 100 2400 200]);
set(gcf,'PaperPositionMode','auto');
ylim([-2 10])

set(gca,'XTickLabel',[]);

%%
figure(5)
ppp2=violin(Index_Norm_iS,'facecolor','red','edgecolor','k',...
    'bw',0.4,'mc',[],'medc','k','facealpha','0.5')
plot(0:25,4*ones(1,size(noise_level,1)+2),'--','linewidth',2)
rectangle('Position',[0 3 25 2],'EdgeColor','none','FaceColor',[0 1 0 0.4])
chi=get(gca, 'Children')
set(gca, 'Children',flipud(chi))
legend('off')
grid on
legend('off')
set(gca,'FontSize',20);
set(gcf,'Position',[100 100 2400 200]);
set(gcf,'PaperPositionMode','auto');
ylim([-2 10])
set(gca,'XTickLabel',[]);

%%
% We will plot the index norm difference of each result. This will tell us
% how many terms got wrongly selected for each method.
figure(6)
hold on
rectangle('Position',[0 -1 25 2],'EdgeColor','none','FaceColor',[0 1 0 0.4])
plot(0:25,0*ones(1,size(noise_level,1)+2),'--','linewidth',4)
violin(Index_Norm_Diff_DL,'facecolor','blue','edgecolor','k','facealpha',1,...
    'bw',0.4,'mc',[],'medc','k')
% chi=get(gca, 'Children')
% set(gca, 'Children',flipud(chi))
grid on
legend('off')
set(gca,'FontSize',34);
%set(gcf,'PaperPositionMode','auto');
ylim([-2 10])
set(gca,'XTickLabel',[]);
set(gcf,'Position',[100 100 3400 400]);
%%
figure(7)
hold on
rectangle('Position',[0 -1 25 2],'EdgeColor','none','FaceColor',[0 1 0 0.4])
plot(0:25,0*ones(1,size(noise_level,1)+2),'--','linewidth',4)
violin(Index_Norm_Diff_iS,'facecolor','red','edgecolor','k','facealpha',1,...
    'bw',0.4,'mc',[],'medc','k')
%chi=get(gca, 'Children')
%set(gca, 'Children',flipud(chi))
grid on
legend('off')
set(gca,'FontSize',34);
set(gcf,'Position',[100 100 2400 400]);
%set(gcf,'PaperPositionMode','auto');
ylim([-2 10])
noise_level(1)=0;

for i=1:length(noise_level)
    Xlabel_tick{1,i}=num2str(noise_level(i));
end

set(gca,'XTick',1:24)
set(gca,'XTickLabel',Xlabel_tick);

%%
% We also interested in which terms got wrongly selected most in SINDy-PI
figure(8)
subplot(2,1,1)
imagesc(Total_Selection_LHS1/length(mat_DL)*100)
colorbar

subplot(2,1,2)
imagesc(Total_Selection_LHS2/length(mat_DL)*100)
colorbar










