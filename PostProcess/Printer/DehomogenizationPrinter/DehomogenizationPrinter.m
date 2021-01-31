classdef untitled3
fPath = '/home/alex/git-repos/Swan/Output';
fName = 'ExperimentingPlot';
%fName = 'LatticeExperimentInputCantileverSymmetricMeshSuperEllipsePDE';

fNameMsh = [fName,'4.flavia.msh'];
s.filePath = fNameMsh;
wM = WrapperMshFile(s);
wM.read();
dM = wM.getDataBase();

s.nElem  = size(dM.connec,1);
s.nNodes = size(dM.coord,1);
fNameRes = [fName,'4.flavia.res'];
s.filePath = fullfile(fPath,fName,fNameRes);
wR = WrapperResFile(s);
wR.read();
dR = wR.getDataBase();

s.connec = dM.connec;
s.coord  = dM.coord;
m = Mesh(s);
s.mesh   = m;
s.fValues = dR.Alpha(:,1);
f = PieceWiseConstantFunction(s);
alpha1 = f.projectToLinearNodalFunction();

s.fValues = dR.Alpha(:,2);
f = PieceWiseConstantFunction(s);
alpha2 = f.projectToLinearNodalFunction();


fOutName = fullfile(fPath,fName,[fName,'Dehomog.txt']);
s.fileName = fOutName;
s.values(:,1) = dR.M1;
s.values(:,2) = dR.M2;
s.values(:,3) = alpha1;
s.values(:,4) = alpha2;
s.values(:,5) = dR.SuperEllipseExponent;
p = DehomogenisationOutputPrinter(s);
p.print();

fOutName = fullfile(fPath,fName,[fName,'DehomogMesh.txt']);
s.fileName = fOutName;
s.mesh     = m;
p = DehomogenisationMeshPrinter(s);
p.print();

end

