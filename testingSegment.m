% Testing segment

clear;
clc;

mesh        = createMesh();
fun         = createReferenceField(mesh);
nLFilter    = createNonLinearFilter(mesh);
fNL         = nLFilter.compute(fun,2);

chi = fun.project('P1');
chi.plot
fNL.plot
gradChi = Grad(chi);
gradRho = Grad(fNL);
plotGradient(gradChi,mesh);
plotGradient(gradRho,mesh);


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

function nLFilter = createNonLinearFilter(m)
s.mesh   = m;
s.theta  = 0;
s.alpha  = 1*m.computeMeanCellSize();
s.beta   = 10*m.computeMeanCellSize();
s.type   = 'Segment';
nLFilter = NonLinearFilter.create(s);
end

function plotGradient(gradVar,mesh)
gradRho= gradVar.project('P1',mesh);
x = mesh.coord(:,1);
y = mesh.coord(:,2);
t  =gradRho.fValues;
ct = (t(:,1));
st = (t(:,2));


n = 4;  % Modify this value to control density
x = x(1:n:end);
y = y(1:n:end);
ct = ct(1:n:end);
st = st(1:n:end);


figure;
quiver(x, y, ct, st, 'AutoScale', 'on', 'LineWidth', 1.5);  % Increase LineWidth for thicker arrows

axis equal;  % Keep aspect ratio equal
%    grid on;
box on;      % Adds a box around the plot
xlim([min(x), max(x)]);
ylim([min(y), max(y)]);
end