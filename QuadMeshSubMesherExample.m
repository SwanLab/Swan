function QuadMeshSubMesherExample

coord = [0 0;1 0;1 1;0 1;2 0;2 1;0 2;1 2;2 2];
connec = [1 2 3 4; 2 5 6 3; 4 3 8 7; 3 6 9 8];
ls = [-0.05 0.2 -0.5 0.1 0.1 -1 1 -0.2 -0.5]';

[connecFull,connecCut,cutElemsQ] = computeConnecCutAndFull(ls,connec);
mInterior = computeMesh(connecFull,coord);

mCut = computeMesh(connecCut,coord);

lastNode = max(connec(:));
s.mesh = mCut;
s.lastNode = lastNode;
thisSubMesher = SubMesher(s);
subMesh   = thisSubMesher.computeSubMesh();
subCutMesh = subMesh;


ls = computeLevelSetInSubCutMesh(thisSubMesher,ls);
connec = subCutMesh.connec;
coord  = subCutMesh.coord;


[connecSubCellFull,connecSubCellCut,cutElems] = computeConnecCutAndFull(ls,connec);



s.coord  = coord;
s.connec = connec;
backgroundCutMesh = Mesh().create(s);
backgroundCutMesh.computeEdges();

cutMeshInteriorTriangle = computeCutMeshInterior(coord,connecSubCellCut,cutElems,ls);
c = cutMeshInteriorTriangle.cellContainingSubcell;

nCut = size(connecCut,1);
nnode = size(connecCut,2);

cellCt = repmat((1:nnode)',1,nCut);
cellCt = cellCt(:);

cellofSubStriangle = repmat(cutElemsQ',nnode,1); 
cellofSubStriangle = cellofSubStriangle(:);

cellContainingSubcell = cellofSubStriangle(c);

xIso = cutMeshInteriorTriangle.xCoordsIso;

t = cellCt(c);

xCoord = thisSubMesher.projectToIsoSubMesh(t,xIso);



connecCutInterior = cutMeshInteriorTriangle.connec;
cutMesh.connec = [connecSubCellFull;connecCutInterior];
cutMesh.coord  = cutMeshInteriorTriangle.coord;
cutMesh.cellContainingSubcell = cellContainingSubcell;
cutMesh.xCoordsIso = xCoord;

mCutInterior = computeMesh(cutMesh.connec,cutMesh.coord);
figure
mInterior.plot();
hold on
mCutInterior.plot();



d = load('cutMeshProvisionalQuad');

n1 = norm(d.connec(:) - cutMesh.connec(:));
n2 = norm(d.xCoordsIso(:) - cutMesh.xCoordsIso(:));
n3 = norm(d.cellContainingSubcell - cutMesh.cellContainingSubcell);
n4 = norm(d.coord(:) - cutMesh.coord(:));

error = n1 + n2 + n3 + n4




end


function cMeshInt = computeCutMeshInterior(coord,connec,cutElems,ls)
s.coord  = coord;
s.connec = connec;
backgroundCutMesh = Mesh().create(s);
backgroundCutMesh.computeEdges();

s.backgroundMesh = backgroundCutMesh;
s.cutElems = cutElems;
s.levelSet = ls;
cMeshInt = CutMeshComputerProvisional(s);

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

function ls = computeLevelSetInSubCutMesh(subMesher,ls)
levelSet = subMesher.projectToSubMesh(ls);
ls = [ls;levelSet];
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