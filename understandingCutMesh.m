function understandingCutMesh
coord = [0 0;1 0;1 1;0 1;2 0;2 1;0 2;1 2;2 2; 0.5 0.5; 1.5 0.5; 0.5 1.5; 1.5 1.5];
connec = [1 2 10; 2 3 10; 10 3 4; 10 4 1; 2 11 3; 2 5 11; 5 6 11; 11 6 3; 3 8 12; 4 3 12; 12 8 7; 12 7 4; 3 6 13; 6 9 13; 13 9 8; 3 13 8];
s.coord = coord;
s.connec = connec;
m = Mesh_Total(s);

ls = [-0.05 0.2 -0.5 0.1 0.1 -1 1 -0.2 -0.5 -0.05 -0.05 0.05 -0.5]';

% inter = Interpolation.create(m,'LINEAR');
% s.unfittedType = 'INTERIOR';
% s.meshBackground = m;
% s.interpolationBackground = inter;
% s.includeBoxContour = true;
% uMesh = UnfittedMesh(s);
% uMesh.compute(ls);
% uMesh.plot()

lsV(:,1) = ls(connec(:,1));
lsV(:,2) = ls(connec(:,2));
lsV(:,3) = ls(connec(:,3));
isFull = all(lsV<0,2);
isEmpty = all(lsV>0,2);
isCut = ~isFull & ~isEmpty;
connecFull = connec(isFull,:);
connecCut = connec(isCut,:);

s.coord = coord;
s.connec = connec;
s.connec = connecCut;

s.levelSet = ls;
cutMesh = CutMeshComputerProvisional(s);
connecCutInterior = cutMesh.connecCutInt;
coordT = cutMesh.coord;



s.connec = [connecFull;connecCutInterior];
s.coord = coordT;
figure
m = Mesh().create(s);
m.plot();
end




