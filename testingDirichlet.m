%% Testing subfunctions
% Goal: generate subfunctions of functions (eg. left boundary)
clear; close all;

% Create P1 Function
mesh = UnitQuadMesh(10,10);
sAF.fHandle = @(x) x(2,:,:);
sAF.ndimf   = 1;
sAF.mesh    = mesh;
xFun = AnalyticalFunction(sAF);

p1 = xFun.project('P1');

% Boundary submesh
isMiddle = @(coor) (abs(coor(:,2) >= 0.4 & abs(coor(:,2)) <= 0.6) );
[mBound, l2gBound] = mesh.createSingleBoundaryMesh();
validNodes = find(isMiddle(mBound.coord));
validElems = find(sum(ismember(mBound.connec, validNodes),2) == 2); % == 2 because line
coord_valid  = mBound.coord(validNodes, :);
connec_valid = mBound.connec(validElems,:);

newNodes = (1:size(coord_valid,1))';
nElemNew = size(connec_valid, 1);
newConnec = [(1:nElemNew-1)', (2:nElemNew)']; % we need to check again if the connectivities are valid!

s.connec = newConnec;
s.coord = coord_valid;
s.kFace = -1;

m = Mesh(s);
m.plot

l2gSub(newNodes(:)) = l2gBound(validNodes);

% Create SubdomainMesh class: mesh, parent mesh + local 2 global
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
