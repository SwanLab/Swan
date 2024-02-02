%% Testing other functions
% Create a Mesh FEM results
clear; close all;

% Create P1 Function
mesh = UnitQuadMesh(10,10);
sAF.fHandle = @(x) x(2,:,:);
sAF.ndimf   = 1;
sAF.mesh    = mesh;
xFun = AnalyticalFunction(sAF);

p1 = xFun.project('P1');


isMiddle = @(coor) (abs(coor(:,2) >= 0.4 & abs(coor(:,2)) <= 0.6) );

% % Boundary mesh
% x = mesh.coord(:,1);
% y = mesh.coord(:,2);
% 
% k = boundary(x,y);
% k = k(1:end-1);
% originalNodes = k;
% newNodes = (1:length(k))';
% boundaryCoords = [x(k), y(k)];
% boundaryConnec = [newNodes, circshift(newNodes,-1)];
% 
% % hmmm
% aa.connec = boundaryConnec;
% aa.coord = boundaryCoords;
% aa.kFace = -1;
% 
% bM = Mesh(aa);
% local2global(newNodes(:)) = originalNodes(:);

