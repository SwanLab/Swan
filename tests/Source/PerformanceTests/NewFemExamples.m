clc; clear;

s.dim    = '2D';
s.length    = 1;
s.height = 0.1;
test = PerformanceTest(s);
sol = test.compute(0.05);
%% 
clc; clear; close all;

% 2D
% index = 1;
% for i = 0.01:0.005:0.5
%     gg = @() example2D(i);
%     temps(index) = timeit(gg);
%     connecs(index) = size(example2D(i),1);
%     index = index+1;
% end
% 
% plot(0.01:0.005:0.5, temps)
% plot(0.01:0.005:0.5, connecs)

% 3D
for i = 2:1:10
    gg = @() example3D(i);
    temps(i) = timeit(gg);
    connecs(i) = size(example3D(i),1);
end

plot(1:1:10, temps)
plot(1:1:10, connecs)

% example3D()

%% 2D Example
% Cantilever beam
function connec = example2D(mstep)
% clc; clear; close all;

% Coordinates
x = 0:  mstep  :1;
y = 0: mstep/2 :0.25;

% Mesh II

s.dim    = '2D';
s.len    = 1;
s.height = 0.1;
beam = CantileverBeam(s);
mesh = beam.create(5,5);
coords = mesh.coord;
connec = mesh.connec;
npnod = size(coords,1);

% Boundary conditions
dirichNodes = find(coords(:,1) == 0 );
ndirich = size(dirichNodes,1);
dirichlet = [dirichNodes,   ones(ndirich,1), zeros(ndirich,1);
             dirichNodes, 2*ones(ndirich,1), zeros(ndirich,1)];
neumann   = [npnod, 2, 0.0001];

% FEM parameters
p.dim = '2D';
p.type = 'ELASTIC';
p.scale = 'MACRO';
p.mesh  = mesh;
p.dirichlet = dirichlet;
p.pointload = neumann;

% Solution
fem = NewFEM.create(p);
fem.solve();
% fem.plot();
end

%% 3D Example
% Cantilever beam
function connec  = example3D(mstep)
% clc; clear; close all;

% Coordinates
% x = 0:  mstep  :1;
% y = 0: mstep/2 :0.25;
% z = 0: mstep/2 :0.25;

% Mesh II
s.dim    = '3D';
s.len    = 1;
s.height = 0.1;
beam = CantileverBeam(s);
mesh = beam.create(5,5);
coords = mesh.coord;
connec = mesh.connec;
npnod = size(coords,1);

% Boundary conditions
dirichNodes = find(coords(:,1) == 0 );
ndirich = size(dirichNodes,1);
dirichlet = [dirichNodes,   ones(ndirich,1), zeros(ndirich,1);
             dirichNodes, 2*ones(ndirich,1), zeros(ndirich,1);
             dirichNodes, 3*ones(ndirich,1), zeros(ndirich,1)];
neumann   = [npnod, 2, 0.0001];

% FEM parameters
p.dim = '3D';
p.type = 'ELASTIC';
p.scale = 'MACRO';
p.mesh = mesh;
p.dirichlet = dirichlet;
p.pointload = neumann;

% Solution
fem = NewFEM.create(p);
fem.solve();
% fem.plot();

end
%% Visualize

%             Tn = connec;
%             x  = coords(:,1);
%             y  = coords(:,2);
%             figure()
%             hold on
%             colormap jet;
%             plot(x(Tn)',y(Tn)','--','linewidth',0.5);
