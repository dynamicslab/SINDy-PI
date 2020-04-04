%% This file is the function file of an example implicit PDE file. We manully code up the best model discovered by the DL-SINDy
% Coded By: K
% Last Updated: 2019/06/17
%% 
function rhs=Ex_ODE_1_Best_Model(t,u0t,k,x)
% Back to time domain
u0=ifft(u0t); 

% ux: In time domain
u0x=ifft(i*k.*u0t);

% uxx: In time domain
u0xx=ifft((i*k).^2.*u0t);

% uxx: In time domain
u0xxx=ifft((i*k).^3.*u0t);

% ux: In Fourier domain
u0x_f=((i*k).^1).*u0t;

% uxx: In Fourier domain
u0xx_f=((i*k).^2).*u0t;

% uxx/(u): In time domain
gx=(5.972e17*u0 - 5.233e15*u0x + 8.452e18*u0xx - 1.387e16*u0xxx + 2.388e15.*u0.*u0x + 1.159e18.*u0.*u0xx - 1.396e17)./(1.723e18*u0 + 2.992e18) ;

% uxx/(u+10): In Fourier domain
Gx_f=fft(gx);

% FFT of RHS
rhs =Gx_f;







