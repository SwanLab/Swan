
clear;

%% INPUT
folder   = 'TFGAlvaro/Results/';
fileName = 'Bicycle_01';
gidName  = 'Bicycle';
testName = 'testAlvaro'; % important to use the same mesh as here !!

%% CODE
refLine = '            <DataArray Name="fValues" NumberOfComponents="3" format="ascii" type="Float64">';
fileID = fopen([folder,fileName,'.vtu']);
tline  = fgetl(fileID);
while ischar(tline)
    switch tline
        case refLine
            tline = fgetl(fileID);
            break;
        otherwise
            tline = fgetl(fileID);
    end
end

a.fileName = gidName;
s = FemDataContainer(a);
mesh = s.mesh;
elParams.mesh = mesh;
elParams.material = s.material;
elParams.scale = 'MACRO';
elParams.bc = s.bc;
physicalProblem = ElasticProblem(elParams);


fVal = zeros(mesh.nnodes,1);
for i = 1:mesh.nnodes
    fk      = find(tline(2:end)==' ')+1;
    fVal(i) = str2double(tline(1:fk(1)-1));
    tline   = fgetl(fileID);
end
fclose(fileID);
sFun.fValues = fVal;
sFun.mesh    = mesh;
dv.fun  = P1Function(sFun);
dv.type = 'LevelSet';
dv.mesh = mesh;
desVar  = LevelSet(dv);

set = Settings(testName);
set.warningHoleBC = false;
set.printIncrementalIter = false;
set.printChangingFilter = false;
set.printing = false;
translator = SettingsTranslator();
translator.translate(set);
set  = SettingsTopOptProblem(testName);
hSet = set.homogenizedVarComputerSettings;
homogenizedVarComputer = HomogenizedVarComputer.create(hSet);

chi           = desVar.getCharacteristicFunction();
ss.filterType = 'LUMP';
ss.mesh       = mesh;
ss.trial      = P1Function.create(mesh,1);
filter        = Filter.create(ss);
fP1Lump       = filter.compute(chi,'QUADRATICMASS');
q             = physicalProblem.getQuadrature();
xP0           = squeeze(fP1Lump.evaluate(q.posgp));
xf            = cell(2,1);
xf{1}         = reshape(xP0',[mesh.nelem,q.ngaus]);
xf{2}         = [ones(1,mesh.nelem);zeros(2,mesh.nelem)];
homogenizedVarComputer.computeCtensor(xf);
%physicalProblem.setC(homogenizedVarComputer.C);
physicalProblem.solve();

%% OUTPUT
physicalProblem.uFun.print([folder,fileName,'Displacements_test'],'Paraview');