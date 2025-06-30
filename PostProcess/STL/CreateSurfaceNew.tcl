proc CreateSurfaceNew {arg1 arg2 arg3 arg4 arg5} {
    set input $arg1
    set output $arg2
    set meshName $arg3
    set gidProjectName $arg4
    set gidBasPath $arg5
    GiD_Process Mescape Postprocess
    GiD_Process Mescape Files Read -rebuildIndex:0 $input
    GiD_Process Mescape Results Options ExtractBoundaries escape
    GiD_Process Mescape Results IsoSurfaces DisplayStyle ContourFillColor escape
    GiD_Process Mescape Results IsoSurfaces Exact fValues 1 0.0
    GiD_Process Mescape 'Redraw
    GiD_Process Mescape Results IsoSurfaces TurnIntoCuts escape
    GiD_Process Mescape DoCut ConvertToSets ALL_CUTSETS escape
    foreach set_name {"Set0 IsoSurf x = 0 WORKPIECE" "WORKPIECE boundary"} {
    GidUtils::CreateMeshFromSet $set_name Layer0
    }
    GiD_Process Mescape files saveall allmeshessets $meshName
    GiD_Process Mescape Preprocess Yes
    GiD_Process Mescape Files MeshRead "$meshName.msh"
    GiD_Process Mescape Files WriteForBAS $gidBasPath $output
    GiD_Process Mescape Files DxfRead -collapse:1 -ignorelayer:1 -tolerance:0.0001 -autotol:1 -- $output
    GiD_Process Mescape Geometry Create NurbsSurface 1:end escape
    GiD_Process Mescape Files Save -alsoresults:1 -- geoversion:current $gidProjectName
    GiD_Process Mescape Geometry Delete AllTypes
}
