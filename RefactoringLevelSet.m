clear all;
%close all;

path = "C:\Users\berna\Desktop\TFG\Swan\TopOptEig\OptimalBucklingColumn\Results";
v = load(path + "\Circular hole\C2.mat");
m = load(path + "\mesh.mat");

value = v.val;

mesh = m.mesh;

% h = 3*ones(500,1);
% b = 2*ones(500,1);
% eh = 0.5*ones(500,1);
% eb = ones(500,1);
% value = [h;b;eh;eb];
s.designVariableValue = value; %valueHole(1:1000)/valueFull(1:500)
s.coord = mesh.coord;
s.type = 'holedCircle'; %'cylinderBuckling'/'holedCircle'/'rectangularColumn'/'rectangularHoleColumn'
plt = Plot3DBucklingColumn(s);
plt.compute();

 