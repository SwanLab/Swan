function QuadBoundaryMeshExperiment
close all
coord = [0 0;1 0;1 1;0 1;2 0;2 1;0 2;1 2;2 2];
connec = [1 2 3 4; 2 5 6 3; 4 3 8 7; 3 6 9 8];
ls = [-0.05 0.2 -0.5 0.1 0.1 -1 1 -0.2 -0.5]';
%ls = [-0.05 0.2 -0.5 -0.1 0.1 -1 1 -0.2 -0.5]';
nameCase = 'cutMeshProvisionalQuad';
s.coord = coord;
s.connec = connec;
m = Mesh().create(s);

figure()
s.meshBackground = m;
s.unfittedType = 'INTERIOR';
uMesh = UnfittedMesh(s);
uMesh.compute(ls);
uMesh.plot();
bMeshOld = uMesh.boundaryCutMesh;

[connecFull,connecCut,cutElems] = computeConnecCutAndFull(ls,connec);
mInterior = computeMesh(connecFull,coord);

sM.coord = coord;
sM.connec = connecCut;
backgroundCutMesh = Mesh().create(sM);

s.backgroundMesh = backgroundCutMesh;
s.cutElems = cutElems;
s.levelSet  = ls;
s.lastNode = max(connec(:));
cutMesh = CutMeshProvisionalQuadrilater(s);
cutMesh.compute();

% figure
% mCutInterior = computeMesh(cutMesh.connec,cutMesh.coord);
% figure(10)
% mInterior.plot();
% hold on
% mCutInterior.plot();


%xG =  uMesh.boundaryCutMesh.computeIsoGaussPoints(quad);

quad = Quadrature.set(uMesh.innerCutMesh.geometryType);
quad.computeQuadrature('QUADRATIC2');
m = uMesh.innerCutMesh.cutMeshOfSubCellLocal;
xG2 = m.computeXgauss(quad.posgp);
%xG3 = reshape(xG2,2,[])';

s.coord = m.coord;
conC = connec(uMesh.innerCutMesh.cellContainingSubcell,:);
s.connec = conC;
m = Mesh().create(s);

xC = uMesh.innerCutMesh.cutMeshOfSubCellGlobal.computeXgauss(xG2);
%figure(10)
hold on
plot(xC(1,:),xC(2,:),'*','MarkerSize',10)


bMesh = cutMesh.computeBoundaryMesh();


% 
% figure
% bMeshOld.plot()
% figure%hold on
%bMesh.plot();
quad = Quadrature.set(bMesh.geometryType);
quad.computeQuadrature('QUADRATIC');

[~,m] = cutMesh.obtainXcutIsoBoundary;
xCutB2 = m.computeXgauss(quad.posgp);

xCutB = bMesh.computeXgauss(quad.posgp);
plot(xCutB(1,:),xCutB(2,:),'g*','MarkerSize',10)


patch('vertices',bMesh.coord,'faces',bMesh.connec,...
    'edgecolor',[0 0 1], 'edgealpha',0.5,'edgelighting','flat',...
    'facecolor',[1 0 0],'facelighting','flat','LineWidth',2)
axis('equal');


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
lsElem(:,4) = ls(connec(:,4));
end