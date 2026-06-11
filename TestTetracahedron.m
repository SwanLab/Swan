clc,clear,close all
baseMesh = TetraMesh(1, 1, 1, 20, 20, 20);

gPar.type = 'Octahedron'; %'Octahedron';
gPar.xCoorCenter = 0.5;
gPar.yCoorCenter = 0.5;
gPar.zCoorCenter = 0.5;
gPar.radius = 3/4;
g         = GeometricalFunction(gPar);
phiFun    = g.computeLevelSetFunction(baseMesh);
lsCircle  = phiFun.fValues;
ls = -lsCircle;


sUm.backgroundMesh = baseMesh;
sUm.boundaryMesh   = baseMesh.createBoundaryMesh;
uMesh              = UnfittedMesh(sUm);
uMesh.compute(ls);
uMesh.plot

ls = CharacteristicFunction.create(uMesh);
s.trial = LagrangianFunction.create(baseMesh,1,'P1');
s.mesh = baseMesh;
f = FilterLump(s);
lsf = f.compute(ls,2);

%% GiD file 
file = 'Tetracaedro';
a.fileName = file;
fem = FemDataContainer(a);
nodes = [1:length(fem.mesh.coord)]';
coordNode = [nodes,fem.mesh.coord];

%% Linking dofs 
faceNodes = fem.newBC.microSlave(:,1);
coordFaces = coordNode(faceNodes,:);

planeCoeffs = [1  0  0   1     1    1     1;
               0  1  0   1     1   -1    -1;
               0  0  1   1    -1    1    -1;
               0  0  0  9/4   5/4  5/4   1/4;
               1  1  1  3/4  -1/4 -1/4  -5/4];
figure()
hold on
MS = [];
for i=1:size(planeCoeffs,2)
    n = planeCoeffs(1:3,i);

    fMaster = (abs(coordFaces(:,2:4)*n - planeCoeffs(4,i)) < 1e-2);
    fSlave  = (abs(coordFaces(:,2:4)*n - planeCoeffs(5,i)) < 1e-2);

    coordMaster = coordFaces(fMaster,:);
    coordSlave  = coordFaces(fSlave,:);
    scatter3(coordMaster(:,2),coordMaster(:,3),coordMaster(:,4));
    scatter3(coordSlave(:,2),coordSlave(:,3),coordSlave(:,4));

    distVec = computeDistanceVector(coordMaster,n);
    coordMasterNew = [coordMaster(:,1), coordMaster(:,2:end) + distVec];
    cM = sortrows(coordMasterNew,[2:4]); cS = sortrows(coordSlave,[2:4]);
    max(abs(cM(:,2:4)-cS(:,2:4)),[],'all') % Check position difference between nodes
    MS = [MS; [cM(:,1) cS(:,1)]];
end

%% Clean edge and vertices
% Remove edges and vertices
vals = MS(:);
[uv,~,idx] = unique(vals);
counts = accumarray(idx,1);
repeatedVals = uv(counts == 1);
rowsToKeep = any(ismember(MS,repeatedVals),2);
MSfaces = MS(rowsToKeep,:);
% Obtain edges
edges = uv(counts == 2);
rowsToKeep = any(ismember(MS,edges),2);
MSedges = MS(rowsToKeep,:);
MSedges = removeAndOrder(MSedges);
%Obtain vertices
edges = uv(counts == 3);
rowsToKeep = any(ismember(MS,edges),2);
MSvertices = MS(rowsToKeep,:);
MSvertices = removeAndOrder(MSvertices);

MSFinal = [MSfaces;MSedges;MSvertices]; 

% Plot faces 
idxFaces    = ismember(coordFaces(:,1),MSfaces(:));
idxEdges    = ismember(coordFaces(:,1),unique(MSedges(:)));
idxVertices = ismember(coordFaces(:,1),unique(MSvertices(:)));

figure()
hold on
scatter3(coordFaces(idxFaces,2),coordFaces(idxFaces,3),coordFaces(idxFaces,4),'MarkerEdgeColor','b');
scatter3(coordFaces(idxEdges,2),coordFaces(idxEdges,3),coordFaces(idxEdges,4),'MarkerEdgeColor','g');
scatter3(coordFaces(idxVertices,2),coordFaces(idxVertices,3),coordFaces(idxVertices,4),'MarkerEdgeColor','r');


%% Functions 

function MSnew = removeAndOrder(MS)
    MSnew = [];
    while(~isempty(MS))
        idx2add = find(MS(:,1) == MS(1,1));
        idx2permute = find(MS(:,2) == MS(1,1));
        MSnew = [MSnew; MS(idx2add,:); [MS(idx2permute,2), MS(idx2permute,1)]];

        unique(MS([idx2add;idx2permute],:));
        vals2remove = unique(MS([idx2add;idx2permute],:));

        idx2remove = any(ismember(MS,vals2remove),2);
        MS(idx2remove,:) = [];
    end
end

function [distVec] = computeDistanceVector(coordFace,n)
    centerVec = coordFace(1,2:4) - [0.5,0.5,0.5];
    distVec = -2*(dot(n,centerVec)/norm(n)^2)*n';
end
