proc ExportMSH {arg1 arg2} {
    set input $arg1
    set gidBasPath $arg2
    GiD_Process Mescape Postprocess
    GiD_Process Mescape Files Read $input
    GiD_Process Mescape files saveall allmeshessets "PostProcess/STL/HmmLetMeCook"
}