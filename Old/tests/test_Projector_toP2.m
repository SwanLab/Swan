% Reset Workspace
clear
close all
clc

% Define number of elements
N = [4 8 16 32 64];

% Prepare analytical function
sAF.fHandle = @(x) sin(2*x(1,:,:)*2*pi); % f(x) = sin(x)
%sAF.fHandle = @(x) x(1,:,:).*x(1,:,:).*x(1,:,:)+x(2,:,:).*x(2,:,:).*x(2,:,:);
sAF.ndimf   = 1;

% Variables to store the error
e_p1 = zeros(size(N));
e_p2 = zeros(size(N));
e_p3 = zeros(size(N));

for i = 1:length(N)
    % Create the mesh of N(i) elements
    sAF.mesh = createMesh(N(i));
    
    % Define the analytical function for the selected mesh
    xFun = AnalyticalFunction(sAF);
    
    % Project to P1 and P2
    p1fun = xFun.project('LINEAR');
    p2fun = xFun.project('QUADRATIC');
    p3fun = xFun.project('CUBIC');
    
    r.mesh = sAF.mesh;
    r.quadType = 'ORDER10';
    r.type = 'Error';
    int = Integrator.create(r);
    
    e_p1(i) = int.compute(p1fun,xFun)/p1fun.computeL2norm();
    e_p2(i) = int.compute(p2fun,xFun)/p2fun.computeL2norm();
    e_p3(i) = int.compute(p3fun,xFun)/p3fun.computeL2norm();
    
%     % Define the quadrature to integrate the norm
%     xG = defineGaussPoints(sAF.mesh);
%     
%     % Evaluate the analytical function, and the projections in the Gauss
%     % points
%     f = xFun.evaluate(xG);
%     p1 = p1fun.evaluate(xG);
%     p2 = p2fun.evaluate(xG);
%     p3 = p3fun.evaluate(xG);
%     
%     % Compute the error with the norm of L^2
%     e_p1(i) = sum(sum((f - p1).^2));
%     e_p2(i) = sum(sum((f - p2).^2));
%     e_p3(i) = sum(sum((f - p3).^2));
end

% ----------------------------------------------------------------------- %
% PLOTS
% ----------------------------------------------------------------------- %
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
set(groot,'defaultAxesTickLabelInterpreter','latex');
figure()
loglog(1./N,e_p1,'-x')
hold on
loglog(1./N,e_p2,'-x')
loglog(1./N,e_p3,'-x')
loglog(1./N,1./N.^2,'-x')
loglog(1./N,1./N.^4,'-x')
loglog(1./N,1./N.^6,'-x')
grid on
% legend('P1 error','P2 error','$h^2$','$h^4$','Location','best')
legend('P1 error','P2 error','P3 error','$h^2$','$h^4$','$h^6$','Location','best')
xlabel('h')
ylabel('Error')
% title('P1 and P2 Projections Validation')
title('P1, P2 and P3 Projections Validation')

%% ----------------------------------------------------------------------- %
% FUNCTIONS
% ----------------------------------------------------------------------- %
function m = createMesh(N)
    % Defines a 2D squared mesh of triangles of N nodes by side in a domain
    % between 0 and 1 in both axes
    x1 = linspace(0,1,N);
    x2 = linspace(0,1,N);
    [xv,yv] = meshgrid(x1,x2);
    [F,V] = mesh2tri(xv,yv,zeros(size(xv)),'x');
    sBg.coord  = V(:,1:2);
    sBg.connec = F;
    figure
    m = Mesh.create(sBg);
end

% function xG = defineGaussPoints(mesh)
%     % Returns the gauss points for an integration of order4 for the given
%     % mesh 
%     quad = Quadrature.set(mesh.type);
%     quad.computeQuadrature('ORDER6');
%     xG = quad.posgp;
% end


