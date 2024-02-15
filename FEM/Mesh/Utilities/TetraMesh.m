function m = TetraMesh(length, height, nx, ny)
    s.dim = '3D';
    s.length = length;
    s.height = height;
    cant = CantileverBeamMeshCreator(s);
    m = cant.create(nx, ny);
end