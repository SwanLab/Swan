% SquareIntoEllipse: Microstructures with level set

clear
clc
close all

filenameVec = ["Circle.mat";"HorizontalEllipse.mat";"VerticalEllipse.mat";"DiagonalEllipse.mat"];
for i=1:length(filenameVec)
    filename = char(filenameVec(i));
    plotMicro(filename);
end

function plotMicro(filename)

load(filename,'Mesh','LSval')
bMesh = Mesh.createBoundaryMesh();
s.backgroundMesh = Mesh;
s.boundaryMesh   = bMesh;
unfMesh          = UnfittedMesh(s);
unfMesh.compute(LSval);
unfMesh.createPlotter();
figure
unfMesh.plotStructureInColor('black');

end