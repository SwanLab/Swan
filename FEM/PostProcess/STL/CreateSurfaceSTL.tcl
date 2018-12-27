proc CreateSurfaceSTL {arg1 arg2 arg3 arg4 arg5} {
    set input $arg1
    set output $arg2
    set meshName $arg3
    set gidProjectName $arg4
    set gidBasPath $arg5
    GiD_Process Mescape Postprocess
    GiD_Process Mescape Files Read -rebuildIndex:0 $input
    GiD_Process Mescape Results IsoSurfaces DisplayStyle ContourFillColor
    GiD_Process Mescape Results IsoSurfaces Exact LevelSet 1 0.0
    GiD_Process Mescape 'Redraw
    GiD_Process Mescape Results isosurfaces TurnIntoCuts
    GiD_Process Mescape DoCut TurnIntoMeshes
    GiD_Process Mescape Results Options ExtractBoundaries
    GiD_Process Mescape Utilities Delete SurfaceSets {WORKPIECE 1} Yes
    GiD_Process Mescape Utilities Delete CutSets {IsoSurf LS = 0 WORKPIECE 1} Yes
    GiD_Process Mescape files saveall allmeshessets $meshName
    GiD_Process Mescape Preprocess Yes
    GiD_Process Mescape Files MeshRead "$meshName.msh"
    GiD_Process Mescape Files WriteForBAS $gidBasPath $output
    GiD_Process Mescape Files DxfRead -collapse:1 -ignorelayer:1 -tolerance:0.0001 -autotol:1 -- $output
    GiD_Process Mescape Geometry Create IntMultLines InvertSelection escape Yes
    GiD_Process Mescape Files Save -alsoresults:1 -- geoversion:current $gidProjectName
    GiD_Process Mescape Geometry Delete AllTypes
}
