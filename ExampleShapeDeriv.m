function ExampleShapeDeriv

addpath(genpath(fileparts(mfilename('fullpath'))))
%%%% Trec les propietats de la malla a partir del FemDataContainer %%%%
file = 'CantileverBeam_Triangle_Linear';
a.fileName = file;
s = FemDataContainer(a);

mesh = s.mesh;
mesh.plot()

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


delta = 1e-25;




iter = 1;
tau = 0.01;
step = 1;
while step > delta
    df   = createDerivative(mesh);
    % c(iter) = computeCost(mesh,quad);
    f = createFunction(mesh) ;
    g= project(f,df,mesh);
    newMesh = updateMesh(mesh,g,tau);
    figure(1)
    clf
    newMesh.plot()
    newMesh = remesh(newMesh);
    % axis([0 2 -0.3 1.4])
    % axis([-0.25 2.25 -0.25 1.25])
    axis([0 2 -0.5 1.3])
    figure(2)
    plot(c)
    step = norm(mesh.coord(:) - newMesh.coord(:))
    mesh = newMesh;
    iter   = iter + 1 ;
end

figure(2)
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

fig1 = figure(2);
hold on;
title("\textbf{Cost vs number of iterations}");
plot(c, 'b', 'LineWidth', 1);
xlabel("$Number\ of\ iterations$");
ylabel("$Cost\  of\  the\  function$ ");
grid on;
grid minor;
box on;
set(gcf, 'units', 'points', 'position', [50,50,676/2,420/2]);
hold off;

set(fig1, 'units', 'points');
pos = get(fig1, 'position');
set(fig1, 'PaperPositionMode', 'Auto', 'PaperUnits', 'points', ...
    'PaperSize',[pos(3), pos(4)]);


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
%%elipse%%
 % s.fHandle = @(x) (x(1,:,:)-1).^2/1.5 + (x(2,:,:)-0.5).^2/0.5 -1;
%%cor%%
 % s.fHandle = @(x) ( (x(1,:,:)-1).^2 + (x(2,:,:)-0.3).^2 - 0.5 ).^3 - (x(1,:,:)-1).^2.*(x(2,:,:)-0.3).^3 ;
%%flor%%
 % s.fHandle = @(x) ( (x(1,:,:)-1).^2 + (x(2,:,:)-0.5).^2 ).^3 - 4*(x(1,:,:)-1).^2.*(x(2,:,:)-0.5).^2 ;

%%PDE%%
a = PDEShapeDerivative() ;
u = a.computeTemperature() ;
u_t = a.createTemperatureTarget() ;
s.fHandle = @(x) u - u_t ;


s.ndimf   = 1;
s.mesh    = mesh;
f = AnalyticalFunction(s);
end




function dfC = createDerivative(mesh)
%%elipse%%
  % s.fHandle = @(x) [2*(x(1,:,:)-1)/1.5; 2*(x(2,:,:)-0.5)/0.5];
%%cor%%
  % s.fHandle = @(x) [ 6*( (x(1,:,:)-1).^2 + (x(2,:,:)-0.3).^2 - 0.5 ).^2.*(x(1,:,:)-1) - 2*(x(1,:,:)-1).*(x(2,:,:)-0.3).^3 ;
  %    6*( (x(1,:,:)-1).^2 + (x(2,:,:)-0.3).^2 - 0.5 ).^2.*(x(2,:,:)-0.3) - 3*(x(1,:,:)-1).^2.*(x(2,:,:)-0.3).^2 ] ;
 %%flor%%
  s.fHandle = @(x) [ 6*( (x(1,:,:)-1).^2 + (x(2,:,:)-0.5).^2 ).^2.*(x(1,:,:)-1) - 8*(x(1,:,:)-1).*(x(2,:,:)-0.5).^2 ;
      6*( (x(1,:,:)-1).^2 + (x(2,:,:)-0.5).^2 ).^2.*(x(2,:,:)-0.5) - 8*(x(2,:,:)-0.5).*(x(1,:,:)-1).^2 ] ;
  s.ndimf   = 2;
  s.mesh    = mesh;
  dfC = AnalyticalFunction(s);

% 
% df{1} = @(x)  2*(x(1,:,:)-1)/1.5;
% df{2} = @(x)  2*(x(2,:,:)-0.5)/0.5;
% %df{1} = @(x) x(1,:,:)./abs(x(1,:,:)-1);
% %df{2} = @(x) x(2,:,:)./abs(x(2,:,:)-0.5);
% for i = 1:length(df)
%     s.fHandle = df{i};
%     s.ndimf   = 1;
%     s.mesh    = mesh;
%     dfC{i} = AnalyticalFunction(s);
% end
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