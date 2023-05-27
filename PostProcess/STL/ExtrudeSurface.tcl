proc ExtrudeSurface {arg1} {
    set gidProjectName $arg1
    GiD_Process Mescape PreProcess
    GiD_Process Mescape Files Read $gidProjectName 
    GiD_Process Mescape Utilities Copy Surfaces DoExtrude Surfaces MaintainLayers Translation FNoJoin 0.0,0.0,0.0 FNoJoin 0.0,0.0,0.16 InvertSelection escape
    GiD_Process Mescape Geometry Create volume 1:end escape escape
    GiD_Process Mescape Files Save -alsoresults:1 -- geoversion:current $gidProjectName
}