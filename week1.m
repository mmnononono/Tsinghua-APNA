system = loadcase('case14');
link_branches = [5, 6, 7, 15, 18, 19, 20];  % 题目给的连支
ref_node = 1;

[A, B] = My_makeAB(system, link_branches, ref_node);

disp('A 矩阵:');
disp(full(A));

disp('B 矩阵:');
disp(full(B));

disp('A*B'' 验证结果:');
disp(full(A*B'));
