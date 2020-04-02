function understandingCutMesh
coord = [0 0;1 0;1 1;0 1;2 0;2 1;0 2;1 2;2 2; 0.5 0.5; 1.5 0.5; 0.5 1.5; 1.5 1.5];
connec = [1 2 10; 2 3 10; 10 3 4; 10 4 1; 2 11 3; 2 5 11; 5 6 11; 11 6 3; 3 8 12; 4 3 12; 12 8 7; 12 7 4; 3 6 13; 6 9 13; 13 9 8; 3 13 8];
ls = [-0.05 0.2 -0.5 0.1 0.1 -1 1 -0.2 -0.5 -0.05 -0.05 0.05 -0.5]';

s.coord = coord;
s.connec = connec;

[connecFull,connecCut] = computeConnecCutAndFull(ls,connec);
mInterior    = computeMesh(connecFull,coord);

cutMeshInterior = computeCutMeshInterior(coord,connecCut,ls);

connecCutInterior = cutMeshInterior.connec;
coordCutInterior  = cutMeshInterior.coord;

mCutInterior = computeMesh(connecCutInterior,coordCutInterior);

mInterior.plot();
hold on
mCutInterior.plot();


d = load('cutMeshProvisional');

n1 = norm(d.connec(:) - cutMeshInterior.connec(:));
n2 = norm(d.xCoordIso(:) - cutMeshInterior.xCoordsIso(:));
error = abs(n1 + n2)
end

function m = computeMesh(connec,coord)
s.connec = connec;
s.coord = coord;
m = Mesh().create(s);
end

function cMeshInt = computeCutMeshInterior(coord,connec,ls)
e = computeEdges(connec);
s.cutEdgesParams.nodesInEdges = e.nodesInEdges;
s.cutEdgesParams.levelSet     = ls;
s.cutCoordParams.coord = coord;
s.cutCoordParams.nodesInEdges = e.nodesInEdges;

s.cutEdgesComputerParams.allNodesinElemParams.finalNodeNumber = size(coord,1);
s.cutEdgesComputerParams.allNodesinElemParams.backgroundConnec = connec;
s.cutEdgesComputerParams.allNodesInElemCoordParams.localNodeByEdgeByElem = e.localNodeByEdgeByElem;
s.cutEdgesComputerParams.edgesInElem = e.edgesInElem;
s.cutEdgesComputerParams.nEdgeByElem = e.nEdgeByElem;
s.interiorSubCellsParams.isSubCellInteriorParams.levelSet = ls;

cMeshInt = CutMeshComputerProvisional(s);
end

function [connecFull,connecCut] = computeConnecCutAndFull(ls,connec)
lsInElem = computeLevelSetInElem(ls,connec);
isFull  = all(lsInElem<0,2);
isEmpty = all(lsInElem>0,2);
isCut = ~isFull & ~isEmpty;
connecFull = connec(isFull,:);
connecCut  = connec(isCut,:);
end

function lsElem = computeLevelSetInElem(ls,connec)
lsElem(:,1) = ls(connec(:,1));
lsElem(:,2) = ls(connec(:,2));
lsElem(:,3) = ls(connec(:,3));
end

function e = computeEdges(connec)
s.nodesByElem = connec;
e = EdgesConnectivitiesComputer(s);
e.compute();
end




