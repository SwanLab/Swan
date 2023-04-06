clc;
clear;
close;

x = linspace(0,1,100);
y = linspace(0,1,100);

[xv,yv] = meshgrid(x,y);
s.coord(:,1) = xv(:);
s.coord(:,2) = yv(:);
s.connec = delaunay(s.coord);
m = Mesh(s);
%m.plot()
mesh = m;

% Compute epsilon
L = mesh.computeCharacteristicLength();
nCells = [16 16]; %??
s.epsilon = L./nCells;


% Compute LevelSet
s.type               = 'xd';
s.mesh               = mesh;
s.ndim               = 2;
s.fracRadius = 0.5;
s.widthV= 0.25;
s.widthH= 0.25;
s.creatorSettings    = s;
s.initialCase    = 'Younes';
lSet = LevelSet(s);


% Get unfitted mesh
uM = lSet.getUnfittedMesh();
uM.plot()
