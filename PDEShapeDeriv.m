function PDEShapeDeriv

clear 
clc
addpath(genpath(fileparts(mfilename('fullpath')))) ;

file = 'CantileverBeam_Triangle_Linear';
a.fileName = file;
s = FemDataContainer(a);


% boundaryMesh = BoundaryMesh(s) ;

s.type = 'THERMAL' ;
fem = FEM.create(s);
% fem.solve();





end