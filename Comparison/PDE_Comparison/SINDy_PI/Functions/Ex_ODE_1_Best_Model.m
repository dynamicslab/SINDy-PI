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

% ux: In Fourier domain
u0x_f=((i*k).^1).*u0t;

% uxx: In Fourier domain
u0xx_f=((i*k).^2).*u0t;

% uxx/(u): In time domain
gx=(0.997744*u0 + 9.97524*u0xx + 0.995305.*u0.*u0xx)./(u0+9.97417);

% uxx/(u+10): In Fourier domain
Gx_f=fft(gx);

% FFT of RHS
rhs =Gx_f;







