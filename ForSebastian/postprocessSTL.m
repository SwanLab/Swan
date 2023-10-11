% Postprocess

clear
clc

testName = 'gripping'; % gripping/inverter

load([testName,'.mat']);
close all

levelSet = computation.designVariable;
uMesh    = levelSet.getUnfittedMesh();
% IMcond   = uMesh.createInnerMeshGoodConditioning();
IMcond   = uMesh.createInnerMesh();