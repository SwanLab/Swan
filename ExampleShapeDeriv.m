function ExampleShapeDeriv

file = 'CantileverBeam_Triangle_Linear';
a.fileName = file;
s = FemDataContainer(a);

mesh = s.mesh;
mesh.plot()

fem = FEM.create(s);
fem.solve();
u = fem.uFun;
u.plot
sP.filename = 'Example1';
sP.type = 'GiD';
u.print(sP)


sAF.fHandle = @(x) (x(1,:,:)-1).^2 + (x(2,:,:)-0.5).^2 -0.5;
sAF.ndimf   = 1;
sAF.mesh    = mesh;
xFun = AnalyticalFunction(sAF);


quad = Quadrature.set(s.mesh.type);
quad.computeQuadrature('LINEAR');

fG  = xFun.evaluate(quad.posgp);
dVG = mesh.computeDvolume(quad);
for iF = 1:xFun.ndimf
intF(iF) = sum(squeeze(fG(iF,:,:)).*dVG');
end

sP.mesh = mesh;
sP.connec = mesh.connec;
h = H1Projector_toP1(sP);
g = h.project(xFun)

% Projector to P1
pp1.mesh   = mesh;
pp1.connec = mesh.connec;
projP1 = Projector_toP1(pp1);
p1fun = projP1.project(xFun);
p1fun.plot();

end