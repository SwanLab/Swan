clc;clear;close all;

gm = multicuboid(0.1,0.005,0.005);
E = 210E9;
nu = 0.3;
rho = 7800;
TipVertex = addVertex(gm,'Coordinates',[0.05,0,0.005]);

sModel = createpde('structural','transient-solid');
sModel.Geometry = gm;
msh = generateMesh(sModel);

meshCoord = msh.Nodes';
meshConnec = delaunayn(msh.Nodes');

s.coord    = meshCoord;
s.connec   = meshConnec;
s.kFace    = sModel.Geometry.NumFaces;
meshFine   = Mesh(s);

figure
pdegplot(sModel,'EdgeLabels','on');
title('Beam model')