% Testing non-linear filters

clear;
clc;

mesh        = createMesh();
fun         = createReferenceField(mesh);
filterPDE   = createPDEFilter(mesh);
fPDE        = filterPDE.compute(fun,2);
nLFilter    = createNonLinearFilter(mesh);
fNL         = nLFilter.compute(fun,2);
%errorCircle = 0.5*Integrator.compute((fPDE-fNL).^2,mesh,2);
%disp(errorCircle);





% Functions
function m = createMesh()
x1       = linspace(0,1,100);
x2       = linspace(0,1,100);
[xv,yv]  = meshgrid(x1,x2);
[F,V]    = mesh2tri(xv,yv,zeros(size(xv)),'x');
s.coord  = V(:,1:2);
s.connec = F;
m        = Mesh.create(s);
end

function fun = createReferenceField(m)
s.fHandle = @(x) 1-heaviside((x(1,:,:)-0.5).^2+(x(2,:,:)-0.5).^2-0.3.^2);
s.ndimf   = 1;
s.mesh    = m;
fun       = AnalyticalFunction(s);
end

function filterPDE = createPDEFilter(m)
s.trial        = LagrangianFunction.create(m,1,'P1');
s.filterType   = 'PDE';
s.mesh         = m;
s.boundaryType = 'Neumann';
s.metric       = 'Isotropy';
filterPDE      = Filter.create(s);
end

function nLFilter = createNonLinearFilter(m)
s.trial  = LagrangianFunction.create(m,1,'P1');
s.mesh   = m;
s.type   = 'Circle';
nLFilter = NonLinearFilter.create(s);
end