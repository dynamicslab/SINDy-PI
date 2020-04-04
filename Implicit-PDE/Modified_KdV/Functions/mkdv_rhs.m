%% This file is the function file of the Schwarzian-Kroteweg-de Vires equation
% Coded By: K
% Last Updated: 2019/06/17
%% 
function [rhs,u0t,u0,u0x,u0xx,u0xxx]=mkdv_rhs(t,u0t,k,x,Omega,g0,alpha)
% Back to time domain
u0=ifft(u0t); 

% ux: In time domain
u0x=ifft(i*k.*u0t);

% uxx: In time domain
u0xx=ifft((i*k).^2.*u0t);

% uxxx: In time domain
u0xxx=ifft((i*k).^3.*u0t);

% ux: In Fourier domain
u0x_f=((i*k).^1).*u0t;

% uxx: In Fourier domain
u0xx_f=((i*k).^2).*u0t;

% uxxx: In Fourier domain
u0xxx_f=((i*k).^3).*u0t;

% uxx^2/ux: In time domain
gx = 2*g0./(1+alpha*u0) - 6*u0.*u0x;

% uxx^2/ux: In Fourier domain
Gx_f=fft(gx);

% FFT of RHS
rhs = -u0xxx_f + Gx_f - Omega*u0t;

% Get the derivative in time
u0t=ifft(rhs);







