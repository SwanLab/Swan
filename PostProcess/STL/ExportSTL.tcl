proc ExportSTL {arg1 arg2 arg3} {
    set input $arg1
    set outputSTL $arg2;
    set basFile $arg3
    GiD_Process Mescape PreProcess
    GiD_Process Mescape Files Read $input 
    GiD_Process Mescape Geometry Create NurbsSurface InvertSelection Mescape    
    GiD_Process Mescape Utilities Copy Surfaces DoExtrude Surfaces MaintainLayers Translation FNoJoin 0.0,0.0,0.0 FNoJoin 0.0,0.0,0.16 InvertSelection escape
    GiD_Process Mescape Meshing Generate Yes 0.11 MeshingParametersFrom=Preferences
    GiD_Process Mescape Files WriteForBAS $basFile $outputSTL
    GiD_Process Mescape Quit
}