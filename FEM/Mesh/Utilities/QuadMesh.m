function m = QuadMesh(length, height, nx, ny)
    s.dim = '2D';
    s.length = length;
    s.height = height;
    cant = CantileverBeamMeshCreator(s);
    m = cant.create(nx, ny);
end