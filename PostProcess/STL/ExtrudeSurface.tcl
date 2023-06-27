proc ExtrudeSurface {arg1 arg2} {
    set gidProjectName $arg1
    set height $arg2
    GiD_Process Mescape PreProcess
    GiD_Process Mescape Files Read $gidProjectName 
    GiD_Process Mescape Utilities Copy Surfaces DoExtrude Surfaces MaintainLayers Translation FNoJoin 0.0,0.0,0.0 FNoJoin 0.0,0.0,$height InvertSelection escape
    GiD_Process Mescape Geometry Create volume 1:end escape escape
    GiD_Process Mescape Utilities SwapNormals Surfaces MakeGroupCoherent 1:end escape Yes escape
    GiD_Process Mescape Files Save -alsoresults:1 -- geoversion:current $gidProjectName
}