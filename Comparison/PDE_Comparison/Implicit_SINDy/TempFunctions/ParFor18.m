function f_x = ParFor18(in1)
%PARFOR18
%    F_X = PARFOR18(IN1)

%    This function was generated by the Symbolic Math Toolbox version 8.1.
%    16-Jul-2019 11:21:53

u = in1(:,1);
ux = in1(:,4);
uxxx = in1(:,6);
f_x = (u.*4.174251485746936e-1-u.*ux.*4.282560722695052e-1+u.*uxxx.*1.412308605859289e1-3.35469705192736e-1)./u;
