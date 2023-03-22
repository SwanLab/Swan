x =linspace(0,1,10);
y =linspace(0,1,10);

[xv,yv] = meshgrid(x,y);
s.coord(:,1) = xv(:);
s.coord(:,2) = yv(:);
s.connec = delaunay(s.coord);
m = Mesh(s);
m.plot()
mesh = m;


% LeveSet
epsilon
s.createlll.epsilon = epsilon;
 


uM = ls.getUnfittedMesh();
uM.plot()


% 