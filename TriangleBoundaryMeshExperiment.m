function TriangleBoundaryMeshExperiment
coord = [0 0;1 0;1 1;0 1;2 0;2 1;0 2;1 2;2 2; 0.5 0.5; 1.5 0.5; 0.5 1.5; 1.5 1.5];
connec = [1 2 10; 2 3 10; 10 3 4; 10 4 1; 2 11 3; 2 5 11; 5 6 11; 11 6 3; 3 8 12; 4 3 12; 12 8 7; 12 7 4; 3 6 13; 6 9 13; 13 9 8; 3 13 8];
ls = [-0.05 0.2 -0.5 0.1 0.1 -1 1 -0.2 -0.5 -0.05 -0.05 0.05 -0.5]';
nameCase = 'cutMeshProvisional';

s.coord = coord;
s.connec = connec;
m = Mesh().create(s);

figure()
s.meshBackground = m;
s.unfittedType = 'INTERIOR';
uMesh = UnfittedMesh(s);
uMesh.compute(ls);
bMeshOld = uMesh.boundaryCutMesh;

[connecFull,connecCut,cutElems] = computeConnecCutAndFull(ls,connec);
mInterior = computeMesh(connecFull,coord);

sM.coord = coord;
sM.connec = connecCut;
backgroundCutMesh = Mesh().create(sM);

s.backgroundMesh = backgroundCutMesh;
s.cutElems = cutElems;
s.levelSet = ls;
cutMesh = CutMeshComputerProvisional(s);
cutMesh.compute();

bMesh = cutMesh.computeBoundaryMesh();

bMeshOld.plot()
figure
bMesh.plot();




end

function m = computeMesh(connec,coord)
s.connec = connec;
s.coord = coord;
m = Mesh().create(s);
end

function [connecFull,connecCut,cutElems] = computeConnecCutAndFull(ls,connec)
lsInElem = computeLevelSetInElem(ls,connec);
isFull  = all(lsInElem<0,2);
isEmpty = all(lsInElem>0,2);
isCut = ~isFull & ~isEmpty;
connecFull = connec(isFull,:);
connecCut  = connec(isCut,:);
cutElems = find(isCut);
end

function lsElem = computeLevelSetInElem(ls,connec)
lsElem(:,1) = ls(connec(:,1));
lsElem(:,2) = ls(connec(:,2));
lsElem(:,3) = ls(connec(:,3));
end
