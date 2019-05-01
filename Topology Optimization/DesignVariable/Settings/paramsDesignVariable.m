mesh  = Mesh_GiD('defaultMeshGiD');
value = ones(size(mesh.coord,1),1);
type  = 'Density';
initialCase = 'full';
levelSetCreatorSettings = SettingsLevelSetCreator();