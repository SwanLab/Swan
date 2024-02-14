% Reset Workspace
clear
close all
clc

% Define number of elements
N = 10;

% Prepare analytical function
sAF.fHandle = @(x) sin(x(1,:,:)*(2*pi)); % f(x) = sin(2*pi*x)
sAF.ndimf   = 1;
sAF.mesh = createMesh(N);
    
% Define the analytical function for the selected mesh
xFun = AnalyticalFunction(sAF);

% Project to P1 and P2
p2fun = xFun.project('LINEAR');
p2fun.plot()

function m = createMesh(N)
    % Defines a 2D squared mesh of triangles of N nodes by side in a domain
    % between 0 and 1 in both axes
    x1 = linspace(0,1,N);
    x2 = linspace(0,1,N);
    [xv,yv] = meshgrid(x1,x2);
    [F,V] = mesh2tri(xv,yv,zeros(size(xv)),'x');
    sBg.coord  = V(:,1:2);
    sBg.connec = F;
    
    sBg.coord  = [0,0;1,0;0,1];
    sBg.connec = [1 2 3];
    m = Mesh.create(sBg);
end