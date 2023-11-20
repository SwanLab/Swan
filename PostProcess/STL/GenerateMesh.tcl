proc GenerateMesh {arg1} {
    set gidProjectName $arg1
    GiD_Process Mescape PreProcess
    GiD_Process Mescape Files Read $gidProjectName  
    GiD_Process Mescape Meshing Generate Yes 100 MeshingParametersFrom=Preferences
    GiD_Process Mescape Files Save -alsoresults:1 -- geoversion:current $gidProjectName
    GiD_Process Mescape Files WriteMesh $gidProjectName.msh
}