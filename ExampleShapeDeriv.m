
% function ExampleShapeDeriv
clear 
clc
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


%%%% Defineixo la funció analítica (el funcional F'(Omega)(Theta)) %%%%
% sAF.fHandle = @(x) (x(1,:,:)-1).^2 + (x(2,:,:)-0.5).^2 -0.5;
% sAF.ndimf   = 1;
% sAF.mesh    = mesh;
% xFun = AnalyticalFunction(sAF);

sAF.fHandle = @(x) [2*(x(1,:,:)-1); 2*(x(2,:,:)-0.5)];
sAF.ndimf   = 2;
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
intF = sum(squeeze(fG(:,:,:)).*dVG(:)') ;

%%%% Trobo la funció g que fa que gTheta=F'(Omega)(Theta) %%%%
sP.mesh = mesh;
sP.connec = mesh.connec;
h = H1Projector_toP1(sP);
g = h.project(xFun) ;
% g.plot() ;

newCoord = zeros(length(mesh.coord),2) ;
delta = 0.001;
convergeix = 0 ;
it = 0 ;
coordAsterisco=sAF.mesh.coord ;

while convergeix == 0
    newCoord = coordAsterisco;
    for i = 1:length(mesh.coord)
    newCoord(i,1) = newCoord(i,1) - g.fValues(i,1)*0.001 ;
    newCoord(i,2) = newCoord(i,2) - g.fValues(i,2)*0.001 ;
    end  
    a.connec = mesh.connec ;
    a.coord = newCoord ;
    m = Mesh(a) ;
    t.mesh = m ;
    h = H1Projector_toP1(t);
    g = h.project(xFun) ;

    A = max(max(abs(newCoord-coordAsterisco))) 
    if A < delta
        convergeix = 1 ;
    else 
        coordAsterisco = newCoord ;
        it = it + 1 ;
    end
end

a.connec = mesh.connec ;
a.coord = newCoord ;
m = Mesh(a) ;
m.plot() ;
% xlim([-10,10])
% ylim([-10,10])





% end