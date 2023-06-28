proc ExtrudeSurface {arg1 arg2 arg3} {
    set input $arg1
    set height $arg2
    set output $arg3
    GiD_Process Mescape PreProcess
    GiD_Process Files MeshRead -createLayers:0 -- /home/ton/Github/Swan/PostProcess/STL/HmmLetMeCook.msh
    GiD_Process Mescape Geometry Create SurfMesh SelectMesh 1:end escape
    GiD_Process Mescape Utilities Copy Surfaces DoExtrude Surfaces MaintainLayers Translation FNoJoin 0.0,0.0,0.0 FNoJoin 0.0,0.0,$height InvertSelection escape
    GiD_Process Mescape Geometry Create volume 1:end escape escape
    GiD_Process Mescape Utilities SwapNormals Surfaces MakeGroupCoherent 1:end escape Yes escape
    GiD_Process Mescape Files Save -alsoresults:1 -- geoversion:current $output
}