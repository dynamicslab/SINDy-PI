clc;clear all;close all;
%%
load('Datas/TrainingData.mat')

%%
dt=0.1;
xValue=0:dt:5;

figure(1)
for j=1:7
for i=1:10
    scatter(xValue,xt((i-1)*51+1:i*51,j),50,'filled',...
        'MarkerFaceColor',0.4*[0.1*j 1-0.1*j 0.5],...
        'MarkerEdgeColor',[0 0 0])
    hold on
end
end
box('on')
set(gcf,'Position',[100 100 600 400]);
set(gcf,'PaperPositionMode','auto');
set(gca,'FontSize',18);

