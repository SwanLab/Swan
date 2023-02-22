clear all;
close all;


fC = load("FullCylinderMMA.mat");
hC = load("HoleCylinderMMA.mat");
rC = load("RectangularColumnData.mat");
m = load("MeshBeamOptimizer.mat");

valueFull = fC.val;
valueHole = hC.val;
valueRectangle = rC.val;
mesh = m.meshSample;


s.designVariableValue = valueRectangle(1:1000); %valueHole(1:1000)/valueFull(1:500)
s.coord = mesh.coord;
plt = Plot3DBucklingColumn(s);
plt.compute();


