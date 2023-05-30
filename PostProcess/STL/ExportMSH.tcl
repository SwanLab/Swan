proc ExportMSH {arg1 arg2} {
    set input $arg1
    set gidBasPath $arg2
    GiD_Process Mescape PreProcess
    GiD_Process Mescape Files Read $input
    GiD_Process Mescape Files WriteForBAS $gidBasPath "HmmLetMeCook"
}