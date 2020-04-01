function QuadMeshSubMesherExample

coord = [0 0;1 0;1 1;0 1;2 0;2 1;0 2;1 2;2 2];
connec = [1 2 3 4; 2 5 6 3; 4 3 8 7; 3 6 9 8];
ls = [-0.05 0.2 -0.5 0.1 0.1 -1 1 -0.2 -0.5]';

[connecFull,connecCut] = computeConnecCutAndFull(ls,connec);
mInterior = computeMesh(connecFull,coord);
figure
mInterior.plot();
hold on


mCut = computeMesh(connecCut,coord);

subMesher = SubMesher();
lastNode = max(connec(:));
subCutMesh = computeSubMesh(subMesher,mCut,lastNode);
subCutMesh.plot();


ls = computeLevelSetInSubCutMesh(subMesher,ls);
connec = subCutMesh.connec;
coord  = subCutMesh.coord;


[connecSubCellFull,connecSubCellCut] = computeConnecCutAndFull(ls,connec);

cutMeshInterior = computeCutMeshInterior(coord,connecSubCellCut,ls);

connecCutInterior = cutMeshInterior.connec;
coordCutInterior  = cutMeshInterior.coord;

mCutInterior = computeMesh([connecSubCellFull;connecCutInterior],coordCutInterior);

figure
mInterior.plot();
hold on
mCutInterior.plot();


end


function cMeshInt = computeCutMeshInterior(coord,connec,ls)
s.coord  = coord;
s.connec = connec;
s.cutEdgesParams.edgesComputer = computeEdges(connec);
s.levelSet = ls;
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

function ls = computeLevelSetInSubCutMesh(subMesher,ls)
levelSet = subMesher.projectToSubMesh(ls);
ls = [ls;levelSet];
end

function subMesh = computeSubMesh(subMesher,mCut,lastNode)
figure
s.mesh = mCut;
s.lastNode = lastNode;
subMesh = subMesher.computeSubMesh(s);
end

function mCut = computeMesh(connec,coord)
s.coord = coord;
s.connec = connec;
mCut = Mesh().create(s);
end

function e = computeEdges(connec)
s.nodesByElem = connec;
e = EdgesConnectivitiesComputer(s);
e.compute();
end