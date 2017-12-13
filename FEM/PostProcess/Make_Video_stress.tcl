proc Make_Video_stress {arg1 arg2 arg3 arg4} {
    set postFileList $arg1
    set output_file_name $arg2
    set Field_2_print $arg3 
    set component $arg4 
    GiD_Process Mescape
    GiD_Process Mescape Files ReadMultiple -rebuildIndex:0 $postFileList
    GiD_Process 'Zoom Frame Mescape
    GiD_Process 'AnimationFile Format GIF 
    GiD_Process 'AnimationFile FramesPerStep 5 Mescape 
    GiD_Process 'AnimationFile Start $arg2 Yes Mescape
    GiD_Process  Utilities Variables PostUpdateWindows No Mescape 
    GiD_Process Results ContOptions ColorRamp Tangent Mescape
    
    foreach iStep [GiD_Info post get all_steps {Elastic Problem}] {
    GiD_Process Results AnalysisSel {Elastic Problem} $iStep Mescape
    GiD_Process Results ContourFill $Field_2_print $component Mescape
    GiD_Process DisplayStyle bodybound Mescape
    GiD_Process 'Hardcopy Options ShowLegends No Mescape
    GiD_Process 'Hardcopy Options ShowAxes No Mescape
    GiD_Process 'Hardcopy Options PrintLogo No Mescape
    GiD_Process 'AnimationFile AddStep Mescape 
    }

    GiD_Process 'AnimationFile End
    GiD_Process Utilities Variables PostUpdateWindows Yes Mescape
}