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

h = 3*ones(500,1);
b = 2*ones(500,1);
eh = 0.5*ones(500,1);
eb = ones(500,1);
value = [h;b;eh;eb];
s.designVariableValue = value; %valueHole(1:1000)/valueFull(1:500)
s.coord = mesh.coord;
s.type = 'rectangularHoleColumn'; %'cylinderBuckling'/'holedCircle'/'rectangularColumn'/'rectangularHoleColumn'
plt = Plot3DBucklingColumn(s);
plt.compute();

 