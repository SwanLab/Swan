proc Make_Video_characteristic {arg1 arg2 arg3 arg4 arg5 arg6} {
    set postFileList $arg1
    set output_file_name $arg2
    set Field_2_print $arg3 
    set component $arg4 
    set min_value $arg5 


    GiD_Process Mescape
    GiD_Process Mescape Files ReadMultiple -rebuildIndex:0 $postFileList

    GiD_Process 'Zoom Frame Mescape

    GiD_Process 'Render Smooth
    GiD_Process Results IsoSurfaces DisplayStyle Body Mescape 
    GiD_Process Results IsoSurfaces DisplayStyle Monochrome Mescape
    GiD_Process Results IsoSurfaces DisplayStyle ChangeMonoCol #666666 Mescape 
    GiD_Process Results IsoSurfaces ShowIsolines No Mescape
    GiD_Process Results ContOptions NumberOfColor 1  Mescape
    GiD_Process Results ContOptions SetMaxOptions OutMaxColor Transparent Mescape
    GiD_Process Results ContOptions SetMinOptions OutMinColor Material Mescape
    GiD_Process 'AnimationFile Format GIF 
    GiD_Process 'AnimationFile FramesPerStep 5 Mescape 
    GiD_Process 'AnimationFile Start $arg2 Yes Mescape
    GiD_Process  Utilities Variables PostUpdateWindows No Mescape  

    foreach iStep [GiD_Info post get all_steps {Elastic Problem}] {
    GiD_Process Results AnalysisSel {Elastic Problem} $iStep Mescape 
    GiD_Process Mescape Utilities transformation create yes Mirror FNoJoin 0.0,3.0,0.0 FNoJoin 2.0,3.0,0.0 FNoJoin 0.0,3.0,2.0 escape escape escape escape escape
    GiD_Process 'Zoom Frame Mescape
    GiD_Process 'Rotate Angle $iStep 30 
    GiD_Process Results IsoSurfaces Exact $Field_2_print 1 0 Mescape 
    GiD_Process escape Results ContourFill $Field_2_print $component Mescape 
    GiD_Process Results contoptions setmaxoptions setvalue 0 Mescape
    GiD_Process results contoptions setminoptions setvalue -1e-32 Mescape

    GiD_Process DisplayStyle bodybound Mescape
    GiD_Process 'Hardcopy Options ShowLegends No Mescape
    GiD_Process 'Hardcopy Options ShowAxes No Mescape
    GiD_Process 'Hardcopy Options PrintLogo No Mescape
    GiD_Process 'AnimationFile AddStep Mescape 
    }

    GiD_Process 'AnimationFile End
    GiD_Process 'Hardcopy Options ShowLegends No Mescape
    GiD_Process 'Hardcopy Options ShowAxes No Mescape
    GiD_Process 'Hardcopy Options PrintLogo No Mescape
    GiD_Process 'Rotate Angle 30 30 
    GiD_Process 'Hardcopy PNG $arg6 Mescape
    GiD_Process Utilities Variables PostUpdateWindows Yes Mescape
    GiD_Process Results ContOptions SetMinOptions MinColor Standard Mescape
    GiD_Process Results ContOptions SetMaxOptions MaxColor Standard Mescape
    GiD_Process Results ContOptions SetMinOptions ResetValue Mescape 
    GiD_Process Results ContOptions SetMaxOptions ResetValue Mescape 
    GiD_Process Results ContOptions NumberOfColor 50  Mescape
}

