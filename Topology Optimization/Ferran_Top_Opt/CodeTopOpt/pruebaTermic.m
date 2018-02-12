function  pruebaTermic
%g = [3 4 0 10 10 0 0 0 10 10];
%g = [2 2 2 2;10 10 0 0;10 0 0 10; 0 1 0 0; 1 1 1 0;1 1 0 1;0 0 1 0; 0 0 0 0;0 0 0 0;0 0 0 0];
g = [2 2 2 2;10 10 0 0;10 0 0 10; 0 10 0 0; 10 10 10 0;1 1 1 1;0 0 0 0; 0 0 0 0;0 0 0 0;0 0 0 0];
% polig=[ 2 4  0  0 10 0 10 10 0 10]';
% gd = [polig];
% ns = char('polig')';
% sf = '(polig)';
% [g,bt] = decsg(gd,sf,ns);

% g = [    2.0000   2.0000    1.0000    1.0000    1.0000
%           0        0   -1.0000    0.0000    0.0000
%      1.0000        0    0.0000    1.0000   -1.0000
%           0   1.0000   -0.0000   -1.0000    1.0000
%           0        0   -1.0000         0   -0.0000
%           0        0    1.0000    1.0000    1.0000
%      1.0000   1.0000         0         0         0
%           0        0         0         0         0
%           0        0         0         0         0
%           0        0    1.0000    1.0000    1.0000
% ];
c = 1;a = 0;f = 1;

pdegplot(g,'EdgeLabels','on','SubdomainLabels','on')

model = createpde(1);
geometryFromEdges(model,g);
applyBoundaryCondition(model,'Edge',1:4,'u',0);
generateMesh(model,'Hmax',0.303);
% [K,M,F] = assema(model,c,a,f);
Nodes = model.Mesh.Nodes;
[K,M,F,Q,G,H,R] = assempde(model,c,a,f);
u = assempde(K,M,F,Q,G,H,R);
pdeplot(model,'xydata',u,'zdata',u)



[x,y] = meshgrid(Nodes(1,:),Nodes(2,:));
%     NodesPlot = [x(:)'; y(:)'];
%     z = compute_caracteristic(NodesPl4)ot);
    % zplot = reshape(u,size(x,1),size(x,2));
     tri = delaunay(x,y);
     trisurf(tri,x,y,u)










end