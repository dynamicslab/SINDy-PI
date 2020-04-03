%% This file is the PDE simulation file of Belousov-Zhabotinsky reaction.
% Coded By: K
% Last Updated: 2019/06/24
%%
function [rhs,x_t,z_t,s_t,u_t,x,z,s,u,x_x,z_x,s_x,u_x,x_y,z_y,s_y,u_y,x_xx,z_xx,s_xx,u_xx,x_yy,z_yy,s_yy,u_yy,x_lap,z_lap,s_lap,u_lap]=...
    BZ_Reaction_PDE(t,xzsut,Kx,Kxx,Ky,Kyy,K22,n,N,Dx,Dz,Ds,Du,q,f,ksi,alpha,beta,gama,ksi2,ksi3,phi,NeedDev)

% Calculate u and v terms
xt=reshape((xzsut(1:N)),n,n);
zt=reshape((xzsut((N+1):(2*N))),n,n);
st=reshape((xzsut(2*N+1:3*N)),n,n);
ut=reshape((xzsut((3*N+1):(4*N))),n,n);

x=real(ifft2(xt));
z=real(ifft2(zt));
s=real(ifft2(st));
u=real(ifft2(ut));

% Reaction Terms
xtrhs=reshape((fft2(  (1/ksi)*( f.*z.*(q-x)./(q+x) + x - x.^2 - beta*x + s )  )),N,1);
ztrhs=reshape((fft2(  x - z - alpha*z + gama*u  )),N,1);
strhs=reshape((fft2(  (1/ksi2)*(beta*x - s + phi*u )  )),N,1);
utrhs=reshape((fft2(  (1/ksi3)*(alpha*z - gama*u )  )),N,1);


rhs=[-(Dx/Dx)*K22.*xzsut(1:N)+xtrhs
     -(Dz/Du)*K22.*xzsut(N+1:2*N)+ztrhs
     -(Ds/Du)*K22.*xzsut(2*N+1:3*N)+strhs
     -(Du/Du)*K22.*xzsut(3*N+1:4*N)+utrhs
     ];

 % If you don't need to extract the value, don't run the following code to
 % speed up
if NeedDev==1
    % Get the derivative you want
    x_x=real(ifft2(reshape((1j*Kx).*xzsut(1:N),n,n)));
    z_x=real(ifft2(reshape((1j*Kx).*xzsut(N+1:2*N),n,n)));
    s_x=real(ifft2(reshape((1j*Kx).*xzsut(2*N+1:3*N),n,n)));
    u_x=real(ifft2(reshape((1j*Kx).*xzsut(3*N+1:4*N),n,n)));
    %
    x_y=real(ifft2(reshape((1j*Ky).*xzsut(1:N),n,n)));
    z_y=real(ifft2(reshape((1j*Ky).*xzsut(N+1:2*N),n,n)));
    s_y=real(ifft2(reshape((1j*Ky).*xzsut(2*N+1:3*N),n,n)));
    u_y=real(ifft2(reshape((1j*Ky).*xzsut(3*N+1:4*N),n,n)));
    %
    x_xx=real(ifft2(reshape((1j*Kxx).^2.*xzsut(1:N),n,n)));
    z_xx=real(ifft2(reshape((1j*Kxx).^2.*xzsut(N+1:2*N),n,n)));
    s_xx=real(ifft2(reshape((1j*Kxx).^2.*xzsut(2*N+1:3*N),n,n)));
    u_xx=real(ifft2(reshape((1j*Kxx).^2.*xzsut(3*N+1:4*N),n,n)));
    %
    x_yy=real(ifft2(reshape((1j*Kyy).^2.*xzsut(1:N),n,n)));
    z_yy=real(ifft2(reshape((1j*Kyy).^2.*xzsut(N+1:2*N),n,n)));
    s_yy=real(ifft2(reshape((1j*Kyy).^2.*xzsut(2*N+1:3*N),n,n)));
    u_yy=real(ifft2(reshape((1j*Kyy).^2.*xzsut(3*N+1:4*N),n,n)));
    %
    x_t=real(ifft2(reshape(rhs(1:N),n,n)));
    z_t=real(ifft2(reshape(rhs(N+1:2*N),n,n)));
    s_t=real(ifft2(reshape(rhs(2*N+1:3*N),n,n)));
    u_t=real(ifft2(reshape(rhs(3*N+1:4*N),n,n)));
    %
    x_lap=real(ifft2(reshape(-K22.*xzsut(1:N),n,n)));
    z_lap=real(ifft2(reshape(-K22.*xzsut(N+1:2*N),n,n)));
    s_lap=real(ifft2(reshape(-K22.*xzsut(2*N+1:3*N),n,n)));
    u_lap=real(ifft2(reshape(-K22.*xzsut(3*N+1:4*N),n,n)));
    
    
end


























