function f_x = ParFor19(in1)
%PARFOR19
%    F_X = PARFOR19(IN1)

%    This function was generated by the Symbolic Math Toolbox version 8.1.
%    16-Jul-2019 11:21:53

u = in1(:,1);
ux = in1(:,4);
uxxx = in1(:,6);
f_x = (u.*4.218041615604307e-1-u.*ux.*4.352377068953501e-1+u.*uxxx.*1.436086616420653e1-3.393151414184103e-1)./u;
