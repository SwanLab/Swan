function ExampleShapeDeriv
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
%f = createFunction(mesh);


%%%% Defineixo la quadratura de la malla %%%%
quad = Quadrature.set(s.mesh.type);
% Trec les dades de la quadratura de gauss (ngaus, weigp, posgp)
quad.computeQuadrature('LINEAR');


delta = 0.01;




iter = 1;
tau = 0.01;
step = 1;
while step > delta
    df   = createDerivative(mesh);
    c(iter) = computeCost(mesh,quad);
    f = createFunction(mesh) ;
    g= project(f,df,mesh);
    newMesh = updateMesh(mesh,g,tau);
    figure(1)
    clf
    newMesh.plot()
    newMesh = remesh(newMesh);
    %pause(1)
    axis([0 2 0 1])
    figure(2)
    plot(c)
    step = norm(mesh.coord(:) - newMesh.coord(:));
    mesh = newMesh;
    iter   = iter + 1 ;
end



end

function nMesh = updateMesh(mesh,g,tau)
coord = mesh.coord;
nCoord(:,1) = coord(:,1) - g.fValues(:,1)*tau ;
nCoord(:,2) = coord(:,2) - g.fValues(:,2)*tau ;
s.connec = mesh.connec;
s.coord  = nCoord;
nMesh = Mesh(s);
end

function m = remesh(mesh)
s.coord  = mesh.coord;
s.connec = delaunay(mesh.coord);
m = Mesh(s);
end


function g= project(f,df,mesh)
t.mesh = mesh ;
h = ShapeDerivProjector(t);
g = h.project(f,df);
end


function f = createFunction(mesh)
s.fHandle = @(x) (x(1,:,:)-1).^2 + (x(2,:,:)-0.5).^2 -0.5;
s.ndimf   = 1;
s.mesh    = mesh;
f = AnalyticalFunction(s);
end

function df = createDerivative(mesh)
s.fHandle = @(x) [2*(x(1,:,:)-1); 2*(x(2,:,:)-0.5)];
s.ndimf   = 2;
s.mesh    = mesh;
df = AnalyticalFunction(s);
end

% function g = project1(df,mesh)
% t.mesh = mesh ;
% h = H1Projector_toP1(t);
% g = h.project(df) ;
% end

function c = computeCost(mesh,quad)
f = createFunction(mesh);
fG  = f.evaluate(quad.posgp);
dVG = mesh.computeDvolume(quad);
c = (squeeze(fG)')*dVG';
end

% end