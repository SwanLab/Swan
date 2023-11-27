function m = UnitQuadMesh(nx, ny)
    s.dim = '2D';
    s.length = 1;
    s.height = 2;
    cant = CantileverBeamMeshCreator(s);
    m = cant.create(nx, ny);
end