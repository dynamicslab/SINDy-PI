function f_x = ParFor15(in1)
%PARFOR15
%    F_X = PARFOR15(IN1)

%    This function was generated by the Symbolic Math Toolbox version 8.2.
%    20-Sep-2019 09:35:38

u = in1(:,1);
ux = in1(:,4);
uxxx = in1(:,6);
f_x = ((ux.*-2.895413584232186e19+uxxx.*2.08074189217307e18+u.*ux.*2.276967144755718e19).*-1.483125065871422e-16)./u;
