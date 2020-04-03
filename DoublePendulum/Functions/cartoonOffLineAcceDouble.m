%This function is used to visulize the pendulum movement.
function cartoonOffLineAcceDouble(x,theta1,theta2,L1,L2,HorizontalRange,VerticalRange,displayrate,figNum,JointCol,CartCol,PenCol)
%Define the cart width and hight here.
set(0,'defaultTextInterpreter','latex')
CartWidth=0.3;
CartHeight=CartWidth*0.618;
Radius=0.04;
N=length(x);
dt=0.001;
figIter=1;
figure(figIter);
pin=round(N/figNum);
%Define the display range.
axis([-HorizontalRange HorizontalRange -VerticalRange VerticalRange]);
%Plot zero line
line([-HorizontalRange;HorizontalRange],[-0.5*CartHeight;-0.5*CartHeight],'Color','k','LineWidth',2)
%Define the lables:
pin2=1;
%Update the figure:
pin3=1;
for j=1:1+pin
    if mod(j,displayrate)==0 
        pin3=pin3+1;
    end
end

for m=1:length(x)
    if mod(m,pin)==0
        figIter=figIter+1;
        figure(figIter)
        axis([-HorizontalRange HorizontalRange -VerticalRange VerticalRange]);
        %Plot zero line
        line([-HorizontalRange;HorizontalRange],[-0.5*CartHeight;-0.5*CartHeight],'Color','k','LineWidth',2)
        pin2=0;
        pin3=1;
        for j=m:m+pin
            if mod(j,displayrate)==0
                pin3=pin3+1;
            end
        end
    end
    if mod(m,displayrate)==0 || m==1
        %Clear all previousimage
        %cla
        hold on
        %Define the transparancy
        Tra=0.05+0.95*(pin2/pin3);
        TraCart=0.1+0.2*(pin2/pin3);
        pin2=pin2+1;
        %Define the position of the cart and plot it
        posCart=[x(m,1)-0.5*CartWidth -0.5*CartHeight CartWidth CartHeight];
        rectangle('Position',posCart,'FaceColor',[CartCol TraCart],'Curvature',[0.1,0.1]);
        
        %Plot the base joint
        posBaseBall=[x(m,1)-Radius -Radius 2*Radius 2*Radius];
        rectangle('Position',posBaseBall,'FaceColor',[JointCol Tra],'Curvature',[1,1]);
        
        %Define the position of the pendulum
        xpos=x(m,1)+L1*sin(theta1(m,1))-Radius;
        ypos=L1*cos(theta1(m,1))-Radius;
        %
        x1pos=x(m,1)+L1*sin(theta1(m,1))+L2*sin(theta2(m,1))-Radius;
        y1pos=L1*cos(theta1(m,1))+L2*cos(theta2(m,1))-Radius;
        %Plot the first pendulum arm
        lineX=[x(m,1);x(m,1)+L1*sin(theta1(m,1))];
        lineY=[0;L1*cos(theta1(m,1))];
        line(lineX,lineY,'LineWidth',1.5,'Color',[PenCol Tra]);
        %Plot the second pendulum
        line1X=[x(m,1)+L1*sin(theta1(m,1));x(m,1)+L1*sin(theta1(m,1))+L2*sin(theta2(m,1))];
        line1Y=[L1*cos(theta1(m,1));L1*cos(theta1(m,1))+L2*cos(theta2(m,1))];
        line(line1X,line1Y,'LineWidth',1.5,'Color',[PenCol Tra]);
        
        %Plot joints
        rectangle('Position',posBaseBall,'FaceColor',[JointCol Tra],'Curvature',[1,1]);
        posBall=[xpos ypos 2*Radius 2*Radius];
        rectangle('Position',posBall,'FaceColor',[JointCol Tra],'Curvature',[1,1]);
        pos1Ball=[x1pos y1pos 2*Radius 2*Radius];
        rectangle('Position',pos1Ball,'FaceColor',[JointCol Tra],'Curvature',[1,1]);
        
        
        
        %Make plot propotional
        set(gca,'DataAspectRatio',[1 1 1])
        set(gcf,'Position',[100 100 560*4 420])
        set(gcf,'PaperPositionMode','auto')
        axis off
        %Draw immediately
        drawnow limitrate
    end
end
