clc;
mpc=loadcase('case39');
addpath('../week4'); 
[~,Y]=My_makeY(mpc);
Y=full(Y);

[L,D,U]=My_makeLDU(Y);
E1=Y-L*D*U;
E1_max=max(abs(E1(:)));
E1_Fro=norm(E1,'fro');
fprintf('验证：Y=LDU\n');
fprintf('|Y-LDU|_max=%.3e\n',E1_max);
fprintf('|Y-LDU|_Fro=%.3e\n',E1_Fro);

Z=My_makeZ_fromLDU(L,D,U);
I=eye(size(Y));
E2=Y*Z-I;
E2_max=max(abs(E2(:)));
E2_Fro=norm(E2,'fro');
fprintf('验证：YZ=I\n');
fprintf('|YZ-I|_max=%.3e\n',E2_max);
fprintf('|YZ-I||_Fro=%.3e\n',E2_Fro);
