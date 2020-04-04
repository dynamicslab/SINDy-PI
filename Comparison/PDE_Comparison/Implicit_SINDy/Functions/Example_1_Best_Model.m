%% This file is used to generate the simulation data of an implicit PDE.
% We manully code up the best model to see the performance of it.
% Coded By: K based on the code of Professor Kutz
% Last Updated: 2019/06/17
function [utsol,usol,x,k]=Example_1_Best_Model(dt,T,L,n,PDE_Plot)
%% Define the parameters for the SKdV equation
% Define tspan
tspan=0:dt:T;

% Get the index of each point
x2=linspace(-L/2,L/2,n+1);
x=x2(1:n);

% Get the frequency
k=(2*pi/L)*[0:n/2-1 -n/2:-1].';
k2=fftshift(k);

% Get the initial value
a1=0.5;
u1=(a1*(sech(0.1*(x-15)).^2)).'; 
%
a2=0.5;
u2=(a2*(sech(0.1*(x+15)).^2)).'; 
%
u=u1+u2;

% Transfer the initial point to spectram domain
ut=fft(u);

% Solve for the PDE
opts = odeset('RelTol',1e-12,'AbsTol',1e-13);
[time,utsol]=ode45(@(time,u0)Ex_ODE_1_Best_Model(time,u0,k,x),tspan,ut,opts);

% Back to time domain
for j=1:length(tspan)
    usol(:,j)=ifft( utsol(j,1:n).' );
end

if PDE_Plot==1
    figure()
    hold on
    %
    surf(x,tspan,real(usol'))
    shading interp
    colormap(parula(100))
    %
    Spacer1=10;
    Spacer2=100;
    %
    ss=surf(x(1:Spacer1:size(x,2)),tspan(1:Spacer2:size(tspan,2)),real(usol(1:Spacer1:size(x,2),1:Spacer2:size(tspan,2))'))
    ss.EdgeColor='black'
    %
    set(gca,'FontSize',18);
    set(gcf,'Position',[100 100 600 400]);
    set(gcf,'PaperPositionMode','auto');
    grid on
    
    view(15,65)
end


