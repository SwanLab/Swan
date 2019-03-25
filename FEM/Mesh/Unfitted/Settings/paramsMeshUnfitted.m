unfittedType = 'INTERIOR';
meshBackground = Mesh_GiD('Bridge_Tetrahedra_Coarse');
interpolationBackground = Interpolation.create(meshBackground,'LINEAR');