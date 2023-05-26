% function ExampleShapeDeriv

addpath(genpath(fileparts(mfilename('fullpath')))) ;


%%%% Trec les propietats de la malla a partir del FemDataContainer %%%%
file = 'CantileverBeam_Triangle_Linear';
a.fileName = file;
s = FemDataContainer(a);

mesh = s.mesh;
% mesh.plot()

% fem = FEM.create(s);
% fem.solve();
% u = fem.uFun;
% u.plot
% sP.filename = 'Example1';
% sP.type = 'GiD';
% u.print(sP)

%%%% Defineixo la funció analítica (el funcional F(Omega)) %%%%
sAF.fHandle = @(x) (x(1,:,:)-1).^2 + (x(2,:,:)-0.5).^2 -0.5;
sAF.ndimf   = 1;
sAF.mesh    = mesh;
xFun = AnalyticalFunction(sAF);

%%%% Defineixo la quadratura de la malla %%%%
quad = Quadrature.set(s.mesh.type);
% Trec les dades de la quadratura de gauss (ngaus, weigp, posgp)
quad.computeQuadrature('LINEAR');

% Avaluo la funció analítica a tota la malla, als punts 
fG  = xFun.evaluate(quad.posgp);
% Calculo el diferencial de volum de cada element de la malla
dVG = mesh.computeDvolume(quad);
% Calculo la integral de la funció analítica del funcional F(Omega)
intF = sum(squeeze(fG(:,:,:)).*dVG(:)) ;


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

% end