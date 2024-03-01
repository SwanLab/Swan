proc ExtrudeSurface {arg1 arg2 arg3} {
    set input $arg1
    set height $arg2
    set output $arg3
    GiD_Process Mescape PreProcess
    GiD_Process Files MeshRead -createLayers:0 -- $input
    GiD_Process escape Utilities Collapse mesh Yes

    GiD_Process 'Layers New Layer1 escape
    GiD_Process 'Layers Off Layer0 escape
    GiD_Process 'Layers Off WORKPIECE escape
    GiD_Process 'Layers ToUse Layer1 escape
    GiD_Process Mescape Meshing CreateBoundary Yes
    GiD_Process 'Layers Delete Layer0 Delete WORKPIECE Yes escape
    GiD_Process Mescape Geometry Edit ReConstruction OneLineForEachElement 1:end escape
    GiD_Process Mescape Geometry Create NurbsSurface 1:end escape escape

    GiD_Process Mescape Utilities Copy Surfaces DoExtrude Surfaces MaintainLayers Translation FNoJoin 0.0,0.0,0.0 FNoJoin 0.0,0.0,$height InvertSelection escape
    GiD_Process Mescape Geometry Create volume 1:end escape escape
    GiD_Process Mescape Utilities SwapNormals Surfaces MakeGroupCoherent 1:end escape Yes escape
    GiD_Process Mescape Files Save -alsoresults:1 -- geoversion:current $output
}