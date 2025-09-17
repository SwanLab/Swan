clear;
close all;

%% CODE
refLine = '            <DataArray Name="fValues" NumberOfComponents="1" format="ascii" type="Float64">';
fileID = fopen('NullSLERPResults/TopOpt/MultiLoadBridge/FinalResults_NoOscillations/gJ10_9Loads_fValues.vtu');
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

x1       = linspace(0,10,400);
x2       = linspace(0,2,80);
[xv,yv]  = meshgrid(x1,x2);
[F,V]    = mesh2tri(xv,yv,zeros(size(xv)),'x');
s.coord  = V(:,1:2);
s.connec = F;
mesh = Mesh.create(s);


fVal = zeros(mesh.nnodes,1);
for i = 1:mesh.nnodes
    if i==mesh.nnodes
        fk      = find(tline(2:end)==' ')+1;
        fVal(i) = str2double(tline(1:fk(1)-1));
    else
        fVal(i) = str2double(tline);
    end
    tline   = fgetl(fileID);
end
fclose(fileID);

sUm.backgroundMesh = mesh;
sUm.boundaryMesh   = mesh.createBoundaryMesh();
uMesh = UnfittedMesh(sUm);
uMesh.compute(fVal);
chi = CharacteristicFunction.create(uMesh);

% sF.fValues = fVal;
% sF.mesh = mesh;
% sF.order = 'P1';
% chi = LagrangianFunction(sF);

sFi.filterType = 'LUMP';
sFi.mesh       = mesh;
sFi.trial      = LagrangianFunction.create(mesh,1,'P1');
f            = Filter.create(sFi);
rhoEps = f.compute(chi,3);