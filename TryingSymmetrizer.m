function TryingSymmetrizer
nx = 2;
ny = 3;
x1 = linspace(0,3,nx);
x2 = linspace(0,2,ny);



[xv,yv] = meshgrid(x1,x2);
[F,V] = mesh2tri(xv,yv,zeros(size(xv)),'x');


alpha = 0.3;
ca = cos(alpha);
sa = sin(alpha);
M = [ca -sa; sa ca];

coord = V(:,1:2);
coordN(:,1) = ca*coord(:,1) - sa*coord(:,2);
coordN(:,2) = sa*coord(:,1) + ca*coord(:,2);

close all
s.coord  = coordN;
s.connec = F;
m = Mesh(s);
m.plot()

s.mesh = m;
s.symmetricLine.vector = [ca;sa];
s.symmetricLine.point = [0;0];
mS = MeshSymmetrizer(s);
mS.compute();
m = mS.symmetricMesh;
figure()
m.plot()
end