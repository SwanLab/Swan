%% Mesh refactoring
clc; clear; close all;

% 1D Mesh
a1.coord  = [0.25, 0; 0.35, 0; 0.75, 0; 1, 0];
a1.connec = [1 2; 2 3; 3 4];
a1.kFace  = -1;
mesh1d = Mesh.create(a1);

% 2D Mesh
    % file = 'test2d_triangle';
    % a2.fileName = file;
    % s = FemDataContainer(a2);
    % mesh2d = s.mesh;
mesh2d = UnitTriangleMesh(7,7);
mesh2d = UnitQuadMesh(7,7);

% 3D Mesh
    % file = 'test3d_hexahedra';
    % a3.fileName = file;
    % s = FemDataContainer(a3);
    % mesh3d = s.mesh;
mesh3d = UnitHexaMesh(7,7,7);
mesh3d = UnitTetraMesh(7,7,7);
