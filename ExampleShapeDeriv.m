% function ExampleShapeDeriv

addpath(genpath(fileparts(mfilename('fullpath')))) ;


%%%% Trec les propietats de la malla a partir del FemDataContainer %%%%
file = 'CantileverBeam_Triangle_Linear';
a.fileName = file;
s = FemDataContainer(a);

mesh = s.mesh;
% mesh.plot()

fem = FEM.create(s);
fem.solve();
% u = fem.uFun;
% u.plot
sP.filename = 'Example1';
sP.type = 'GiD';
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

%%%% Trobo la funció g que fa que gTheta=F'(Omega)(Theta) %%%%
sP.mesh = mesh;
sP.connec = mesh.connec;
h = H1Projector_toP1(sP);
g = h.project(xFun) ;
% g.plot() ;

new_coord = zeros(length(mesh.coord),2) ;

for i = 1:length(mesh.coord)
    new_coord(i,1) = sAF.mesh.coord(i,1) - g.fValues(i)*0.01 ;
    new_coord(i,2) = sAF.mesh.coord(i,2) - g.fValues(i)*0.01 ;
end 

mesh2 = mesh ;
mesh2.coord = new_coord ;
mesh2.plot() ;





% end