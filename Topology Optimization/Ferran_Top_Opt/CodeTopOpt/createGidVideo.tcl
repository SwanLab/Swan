proc VideoPrueba {path_Video arg2} {
    #    set pathVideo "/home/aferrerferrer/Desktop/prueba"  #[PasspathVideoCreation]
    set pathVideo "${path_Video}.avi" 
    set postFileList $arg2
#[PassListFiles]
    GiD_Process Mescape Files ReadMultiple -rebuildIndex:0 $postFileList escape escape escape escape escape
    GiD_Process 'Zoom Frame Mescape
    GiD_Process 'AnimationFile Format AVIMSVC
    GiD_Process 'AnimationFile FramesPerStep 5 escape escape escape escape
    GiD_Process 'AnimationFile Start $pathVideo Yes escape escape escape escape escape escape escape escape 
    GiD_Process Utilities Variables PostUpdateWindows no escape escape escape escape Mescape 
    
    foreach iStep [GiD_Info post get all_steps time] {
        GiD_Process Results AnalysisSel time $iStep Mescape Mescape
        GiD_Process Results ContOptions SetMinOptions ResetValue Mescape 
        GiD_Process Results ContOptions SetMaxOptions ResetValue Mescape 
        GiD_Process Results ContourFill GAUSS_YOUNG GAUSS_YOUNG 
        GiD_Process 'AnimationFile AddStep escape escape escape escape 
    }
    GiD_Process Utilities Variables PostUpdateWindows yes escape escape escape escape 
    GiD_Process 'AnimationFile End
}
