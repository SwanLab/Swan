proc ExportSTL {arg1 arg2 arg3} {
    set input $arg1
    set outputSTL $arg2;
    set basFile $arg3
    GiD_Process Mescape PreProcess
    GiD_Process Files MeshRead -createLayers:0 -- $input
    GiD_Process Layers New Layer1 escape
    GiD_Process Layers ToUse Layer0 escape
    GiD_Process Mescape Meshing CreateBoundary Yes
    GiD_Process Layers Delete WORKPIECE Yes escape
    GiD_Process Mescape utilities SwapNormals Select 1:END escape
    GiD_Process Mescape Files WriteForBAS $basFile $outputSTL
}