%%
clc;close all;clear all;
%%
load('Result1.mat')
load('Result2.mat')
load('Result3.mat')

xValue_1=[0;xValue_1;100];
yValue_1=[0;yValue_1;100];

xValue_2=[0;xValue_2;100];
yValue_2=[0;yValue_2;100];

xValue_3=[0;xValue_3];
yValue_3=[0;yValue_3];

LineColor={'b','b','r'};
MKR={'o','^','o'};
%%
figure(1)
hold on

for j=3:-1:1
    if j==3
        qq=3;
    elseif j==2
        qq=2;
    else
        qq=1;
    end
    xValue=eval(strcat('xValue_',num2str(qq)));
    yValue=eval(strcat('yValue_',num2str(qq)));
    
    for i=1:length(yValue)
        if yValue(i)==100
            break
        end
    end
    
    plot([xValue(i) xValue(i)],[0 100],'linewidth',4,'linestyle','-.','color',(LineColor{j}))
    
end

%%
for j=3:-1:1
    if j==3
        qq=3;
    elseif j==2
        qq=2;
    else
        qq=1;
    end
    xValue=eval(strcat('xValue_',num2str(qq)));
    yValue=eval(strcat('yValue_',num2str(qq)));
    
    d_step=0.01;
    p=pchip(xValue,yValue,min(xValue):d_step:max(xValue));
    
    plot(min(xValue):d_step:max(xValue),p,'linewidth',4,'color',(LineColor{j}),'linestyle','-')
    scatter(xValue,yValue,200,MKR{j},'filled',...
        'MarkerFaceColor',(LineColor{j}),'linewidth',2.5,...
        'MarkerEdgeColor',[0 0 0])
    
    
    grid on
    % title('Success Rate vs Data Usage','FontSize',18)
    % xlabel('Data Usage','FontSize', 18)
    % ylabel('Success Rate','FontSize', 18)
    set(gcf,'Position',[100 100 2200 400]);
    set(gcf,'PaperPositionMode','auto');
    set(gca,'FontSize',28);
    
end

%%
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




