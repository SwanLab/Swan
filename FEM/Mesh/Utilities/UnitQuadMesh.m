function m = UnitQuadMesh(nx, ny)
    s.dim = '2D';
    s.length = 1;
    s.height = 1;
    cant = CantileverBeamMeshCreator(s);
    m = cant.create(nx, ny);
end