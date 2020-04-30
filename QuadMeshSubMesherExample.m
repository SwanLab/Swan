function QuadMeshSubMesherExample
close all
coord = [0 0;1 0;1 1;0 1;2 0;2 1;0 2;1 2;2 2];
connec = [1 2 3 4; 2 5 6 3; 4 3 8 7; 3 6 9 8];
ls = [-0.05 0.2 -0.5 0.1 0.1 -1 1 -0.2 -0.5]';
nameCase = 'cutMeshProvisionalQuad';

[connecFull,connecCut,cutElems] = computeConnecCutAndFull(ls,connec);
mInterior = computeMesh(connecFull,coord);

sM.coord  = coord;
sM.connec = connecCut;
backgroundCutMesh = Mesh().create(sM);

s.backgroundMesh = backgroundCutMesh;
s.cutElems = cutElems;
s.levelSet  = ls;
s.lastNode = max(connec(:));
cutMesh = CutMeshProvisionalQuadrilater(s);
cutMesh.compute();

mCutInterior = computeMesh(cutMesh.connec,cutMesh.coord);
figure
mInterior.plot();
hold on
mCutInterior.plot();

d = load(nameCase);

n1 = norm(d.connec(:) - cutMesh.connec(:));
n2 = norm(d.xCoordsIso(:) - cutMesh.xCoordsIso(:));
n3 = norm(d.cellContainingSubcell - cutMesh.cellContainingSubcell);
n4 = norm(d.coord(:) - cutMesh.coord(:));

error = n1 + n2 + n3 + n4

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
lsElem(:,4) = ls(connec(:,4));
end

function mCut = computeMesh(connec,coord)
s.coord = coord;
s.connec = connec;
mCut = Mesh().create(s);
end