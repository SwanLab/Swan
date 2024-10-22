function smallEX
xv = linspace(0,1,2);
yv = linspace(0,1,2);
[X,Y] = meshgrid(xv,yv);
s.coord(:,1) = X(:);
s.coord(:,2) = Y(:);
[F,V] = mesh2tri(X,Y,zeros(size(X)),'f');
s.coord  = V(:,1:2);
s.connec = F;

m = Mesh.create(s);
m.plot

s.mesh = m;
s.fValues(:,1) = m.coord(:,1);
s.fValues(:,2) = -m.coord(:,2);
s.order = 'P1';
p1 = LagrangianFunction(s);
% p1.plot()
p1D = p1.project('P1D');
p1D.plot()
% p11 = p1D.project('P1');
% p11.plot()

mr = m.remesh();
p1Dr = p1D.refine(mr);
p1Dr.plot
p1DrC = p1Dr.project('P1');
p1DrC.plot()
end