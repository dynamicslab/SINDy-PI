%This function is used to visulize the pendulum movement.
function cartoonOffLineAcce(x,theta,L,HorizontalRange,VerticalRange,displayrate)
    %Define the cart width and hight here.
    CartWidth=0.5;
    CartHeight=0.2;
    Radius=0.1;
    %Define the display range.
    axis([-HorizontalRange HorizontalRange -VerticalRange VerticalRange]);
    %Define the lables:
    xlabel('Position')
    title('Anamation')
    %Update the figure:
    for m=1:length(x)
        if mod(m,displayrate)==0
            %Clear all previousimage
            cla
            %Plot zero line
            line([-HorizontalRange;HorizontalRange],[0;0],'Color','k','LineWidth',1.5)
            %Define the position of the cart
            posCart=[x(m,1)-0.5*CartWidth -0.5*CartHeight CartWidth CartHeight];
            rectangle('Position',posCart,'FaceColor',[0.7 .05 .5],'Curvature',[0.1,0.1]);
            %Define the position of the pendulum
            xpos=x(m,1)+L*sin(theta(m,1))-Radius;
            ypos=L*cos(theta(m,1))-Radius;
            posBall=[xpos ypos 2*Radius 2*Radius];
            rectangle('Position',posBall,'FaceColor',[0.1 .02 .3],'Curvature',[1,1]);
            lineX=[x(m,1);x(m,1)+L*sin(theta(m,1))];
            lineY=[0;L*cos(theta(m,1))];
            line(lineX,lineY,'LineWidth',1.5);
            %Make plot propotional
            set(gca,'DataAspectRatio',[1 1 1])
            %Draw immediately
            drawnow limitrate 
        end 
    end
            