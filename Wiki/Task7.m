%% Swan - Wiki
% The student's guide to clean code development
% Task 7: UML of FEM code in Swan repository

% Instructions: run the following code, selecting previously the 'Swan'
% main folder as your current matlab path

file = 'test2d_triangle';
a.fileName = file;
s = FemDataContainer(a);
fem = FEM.create(s);
fem.solve();
fem.print(file);