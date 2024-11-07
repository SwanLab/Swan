
x1 = linspace(0,1,100);
x2 = linspace(0,1,100);
[xv,yv] = meshgrid(x1,x2);
[F,V] = mesh2tri(xv,yv,zeros(size(xv)),'x');
sBg.coord  = V(:,1:2);
sBg.connec = F;
mesh = Mesh.create(sBg);

sAF.fHandle = @(x) [sin(2*pi*x(1,:,:).*x(2,:,:));sin(pi*x(2,:,:).*x(1,:,:))];
sAF.ndimf   = 2; % number of dimensions
sAF.mesh    = mesh;
xFun = AnalyticalFunction(sAF);


p1fun = xFun.project('P1');
%p1fun.plot

t = Voigt(Grad(p1fun)).project('P1',mesh);
%t.plot()
%shading interp