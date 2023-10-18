% Postprocess

clear
clc

% grippingPostPDEFilter, inverterPrePDEFilter
% testName = 'grippingPostPDEFilter'; % gripping/inverter
testName = 'inverterPrePDEFilter'; % gripping/inverter

load([testName,'.mat']);
close all

levelSet = computation.designVariable;
uMesh    = levelSet.getUnfittedMesh();
% IMcond   = uMesh.createInnerMeshGoodConditioning();
IM   = uMesh.createInnerMesh();
EM = IM.computeCanonicalMesh.provideExtrudedMesh(0.05);

