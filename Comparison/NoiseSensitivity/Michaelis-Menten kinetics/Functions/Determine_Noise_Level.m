%% This function will determine the noise level
% Code By: K
% Last Updated: 2019/07/11
%%
function noise_level=Determine_Noise_Level(percent_iter)
% Set the noise level
if percent_iter==1
    noise_level=0;
elseif percent_iter==2
    noise_level=1e-7;
elseif percent_iter==3
    noise_level=5e-7;
elseif percent_iter==4
    noise_level=1e-6;
elseif percent_iter==5
    noise_level=5e-6;
elseif percent_iter==6
    noise_level=1e-5;
elseif percent_iter==7
    noise_level=5e-5;
elseif percent_iter==8
    noise_level=1e-4;
elseif percent_iter==9
    noise_level=5e-4;
elseif percent_iter==10
    noise_level=1e-3;
elseif percent_iter<=14
    noise_level=2e-3*(percent_iter-10);
elseif percent_iter==15
    noise_level=1e-2;
elseif percent_iter<=19
    noise_level=2e-2*(percent_iter-15);
elseif percent_iter==20
    noise_level=1e-1;
elseif percent_iter<=24
    noise_level=1e-1+1e-1*(percent_iter-20);
end

