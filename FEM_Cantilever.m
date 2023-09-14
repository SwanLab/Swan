clc;clear;close all

addpath(genpath(fileparts(mfilename('fullpath'))))

% file = 'Catilever_Mario';
% a.fileName = file;
% s = FemDataContainer(a);

fem = FEM.create(s);
fem.solve();

figure(1)
fem.uFun.plot()
figure(2)
fem.stressFun.plot()
figure(3)
fem.strainFun.plot()