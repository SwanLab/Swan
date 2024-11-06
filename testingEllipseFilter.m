% Testing non-linear filters

clear;
clc;

mesh        = createMesh();
fun         = createReferenceField(mesh);
filterPDE   = createPDEFilter(mesh);
fPDE        = filterPDE.compute(fun,2);
nLFilter    = createNonLinearFilter(mesh);
fNL         = nLFilter.compute(fun,2);
errorEllipse = 0.5*Integrator.compute((fPDE-fNL).^2,mesh,2);
disp(errorEllipse);





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
ss.filterType      = 'PDE';
ss.mesh            =  m;
ss.boundaryType    = 'Neumann';
ss.metric          = 'Anisotropy';
nu                 = 85;    % deg VARIABLE
ss.aniAlphaDeg     = 45;    % alpha FIXED
epsilon            = 1*m.computeMeanCellSize();    % filter radius VARIABLE

ss.CAnisotropic    = [tand(nu), 0; 0, 1/tand(nu)];    % A matrix
ss.trial           = LagrangianFunction.create(m, 1, 'P1');
filterPDE          = Filter.create(ss);
filterPDE.updateEpsilon(epsilon);
end

function nLFilter = createNonLinearFilter(m)
s.trial  = LagrangianFunction.create(m,1,'P1');
s.mesh   = m;
s.type   = 'Ellipse';
theta    = 45;
alpha    = 0.9;
A11 = cosd(theta)^2 + sind(theta)^2 / alpha^2;
A12 = (cosd(theta) * sind(theta) * (alpha^2 - 1)) / alpha^2;
A21 = (cosd(theta) * sind(theta) * (alpha^2 - 1)) / alpha^2;
A22 = sind(theta)^2 + cosd(theta)^2 / alpha^2;
A = [A11, A12; A21, A22];
s.A = inv(A);
s.A = A;

nLFilter = NonLinearFilter.create(s);
end