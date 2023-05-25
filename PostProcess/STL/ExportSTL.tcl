proc ExportSTL {arg1 arg2 arg3} {
    set input $arg1
    set outputSTL $arg2;
    set basFile $arg3
    GiD_Process Mescape PreProcess
    GiD_Process Mescape Files Read $input  
    GiD_Process Mescape Files WriteForBAS $basFile $outputSTL
    GiD_Process Mescape Quit
}