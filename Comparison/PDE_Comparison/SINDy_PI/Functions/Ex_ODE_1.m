%% This file is the function file of an example implicit PDE file
% Coded By: K
% Last Updated: 2019/06/17
%% 
function rhs=Ex_ODE_1(t,u0t,k,x)
% Back to time domain
u0=ifft(u0t); 

% ux: In time domain
u0x=ifft(i*k.*u0t);

% uxx: In time domain
u0xx=ifft((i*k).^2.*u0t);

% ux: In Fourier domain
u0x_f=((i*k).^1).*u0t;

% uxx: In Fourier domain
u0xx_f=((i*k).^2).*u0t;

% uxx/(u): In time domain
gx=u0./(u0+10);

% uxx/(u+10): In Fourier domain
Gx_f=fft(gx);

% FFT of RHS
%rhs =u0x_f+(3/2)*Gx_f;
rhs =u0xx_f+Gx_f;






