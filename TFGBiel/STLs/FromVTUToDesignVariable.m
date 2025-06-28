

clear;
clc;
close all;


file1 = 'TFGBiel/STLs/VTUs/CantileverLevelSet';
mesh1 = TriangleMesh(2,1,150,75);
ls1   = getFValuesFromVTU(file1,mesh1);

s1.backgroundMesh = mesh1;
s1.boundaryMesh   = mesh1.createBoundaryMesh();
uMesh1 = UnfittedMesh(s1);
uMesh1.compute(ls1);






file2 = 'TFGBiel/STLs/VTUs/GripperLevelSet';
gidFile = 'Gripping';
a.fileName = gidFile;
s = FemDataContainer(a);
mesh2 = s.mesh;
ls2 = getFValuesFromVTU(file2,mesh2);

setPositiveLs = @(coor) coor(:,1)>=0.1 & coor(:,1)<=0.13125 & coor(:,2)>=0.45625 & coor(:,2)<=0.54375;
posNodes      = find(setPositiveLs(mesh2.coord));
ls2(posNodes) = 1;

s2.backgroundMesh = mesh2;
s2.boundaryMesh   = mesh2.createBoundaryMesh();
uMesh2 = UnfittedMesh(s2);
uMesh2.compute(ls2);

innerMesh2 = uMesh2.createInnerMesh();
iM3D_2     = innerMesh2.provideExtrudedMesh(0.1);





function fVal = getFValuesFromVTU(file,mesh)
    refLine = '            <DataArray Name="fValues" NumberOfComponents="1" format="ascii" type="Float64">';
    fileID = fopen([file,'.vtu']);
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
    
    fVal = zeros(mesh.nnodes,1);
    for i = 1:mesh.nnodes
        fk      = find(tline(2:end)==' ')+1;
        fVal(i) = str2double(tline(1:fk(1)-1));
        tline   = fgetl(fileID);
    end
    fclose(fileID);
end