function QuadMeshSubMesherExample



coord = [0 0;1 0;1 1;0 1;2 0;2 1;0 2;1 2;2 2];
connec = [1 2 3 4; 2 5 6 3; 4 3 8 7; 3 6 9 8];
%s.coord = coord;
%s.connec = connec;
%m = Mesh().create(s);
%m.plot;

ls = [-0.05 0.2 -0.5 0.1 0.1 -1 1 -0.2 -0.5]';


lsV(:,1) = ls(connec(:,1));
lsV(:,2) = ls(connec(:,2));
lsV(:,3) = ls(connec(:,3));
isFull = all(lsV<0,2);
isEmpty = all(lsV>0,2);
isCut = ~isFull & ~isEmpty;
connecFull = connec(isFull,:);
connecCut  = connec(isCut,:);


s.coord = coord;
s.connec = connecCut;
mCut = Mesh().create(s);
%mCut.plot;


figure
s.mesh = mCut;
s.lastNode = max(connec(:));
subMesher = SubMesher();
subMesh = subMesher.computeSubMesh(s);

subMesh.plot();
levelSet = subMesher.projectToSubMesh(ls);
ls = [ls;levelSet];


lsV2(:,1) = ls(subMesh.connec(:,1));
lsV2(:,2) = ls(subMesh.connec(:,2));
lsV2(:,3) = ls(subMesh.connec(:,3));
isFull = all(lsV2<0,2);
isEmpty = all(lsV2>0,2);
isCut = ~isFull & ~isEmpty;
connecFull2 = subMesh.connec(isFull,:);
connecCut2 = subMesh.connec(isCut,:);

sC.coord  = subMesh.coord;
sC.connec = connecCut2;
sC.cutEdgesParams.edgesComputer = computeEdges(connecCut2);


sC.levelSet = [ls;levelSet];

cutMesh = CutMeshComputerProvisional(sC);
connecCutInterior = cutMesh.connec;
coordT = cutMesh.coord;

s.connec = [connecFull2;connecCutInterior];
s.coord  = coordT;
figure
mCutInterior = Mesh().create(s);
mCutInterior.plot();

hold on

s.connec = connecFull;
s.coord  = coord;
mInterior = Mesh().create(s);
mInterior.plot();


end

function e = computeEdges(connec)
s.nodesByElem = connec;
e = EdgesConnectivitiesComputer(s);
e.compute();
end