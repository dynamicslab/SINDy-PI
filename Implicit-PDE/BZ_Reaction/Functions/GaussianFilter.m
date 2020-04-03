%% This file will generate a Gaussian filter output
% Coded By: K
% Last Updated: 2019/06/25
%%
function weight=GaussianFilter(x,y,wx,wy,fx,fy)
weight=exp(-wx*(x-fx).^2-wy*(y-fy).^2);