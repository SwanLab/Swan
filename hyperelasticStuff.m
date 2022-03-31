filename = 'test_hyperelastic';
s = createFEMparameters(filename);

fem = FEM.create(s);

function s = createFEMparameters(file)
    gidParams = createGiDparameters(file);
    s.dim       = gidParams.pdim;
    s.type      = gidParams.ptype;
    s.scale     = gidParams.scale;
    s.mesh      = gidParams.mesh;
    s.dirichlet = gidParams.dirichlet;
    s.pointload = gidParams.pointload;
end

function gidParams = createGiDparameters(file)
    gidReader = FemInputReader_GiD();
    gidParams = gidReader.read(file);
end