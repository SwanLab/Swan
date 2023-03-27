clc;
clear;
close;

x = linspace(0,1,10);
y = linspace(0,1,10);

[xv,yv] = meshgrid(x,y);
s.coord(:,1) = xv(:);
s.coord(:,2) = yv(:);
s.connec = delaunay(s.coord);
m = Mesh(s);
%m.plot()
mesh = m;

% Compute epsilon
L = mesh.computeCharacteristicLength();
nCells = [10 10]; %??
s.epsilon = L./nCells;


% Compute LevelSet
s.type               = 'Younes';
s.mesh               = mesh;
s.ndim               = 2;
s.creatorSettings    = s;
s.initialCase    = 'Younes';
lSet = LevelSet(s);


% Get unfitted mesh
uM = lSet.getUnfittedMesh();
uM.plot()
