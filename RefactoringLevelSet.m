clear all;
%close all;


fC = load("FullCylinderMMA.mat");
hC = load("HoleCylinderMMA.mat");
rC = load("RectangularColumnData.mat");
rS = load("RectangularColumnSmooth.mat");
m = load("MeshBeamOptimizer.mat");

valueFull = fC.val;
valueHole = hC.val;
valueRectangle = rC.val;
valueRectangleSmooth = rS.val;
mesh = m.meshSample;


s.designVariableValue = valueRectangleSmooth(1:1000); %valueHole(1:1000)/valueFull(1:500)
s.coord = mesh.coord;
plt = Plot3DBucklingColumn(s);
plt.compute();


