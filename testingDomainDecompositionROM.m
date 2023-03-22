% Testing domain decomposition for reduced order modeling

clear
clc


% STEP 1: Full DoF solver

fileName = 'CantileverAxialLoad';

s.testName    = [fileName,'.m'];
s.problemCase = 'cantilever';
s.x1          = 1;
s.y1          = 0.5;
s.N           = 20;
s.M           = 10;
s.P           = 1;
s.DoF         = 2;

CantileverAxial = FEMInputWriter(s);
CantileverAxial.createTest();

a.fileName = fileName;
s = FemDataContainer(a);
fem = FEM.create(s);
fem.solve();
% fem.print(fileName);

u = fem.uFun{1,1};
u.plot();
e = fem.strainFun{1,1};
% e.plot();
sig = fem.stressFun{1,1};
% sig.plot();

% ...


% STEP 2: Mesh decomposition
% ...
