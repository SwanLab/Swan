proc VideoPrueba {def path_Video arg2 arg3} {
    #    set pathVideo "/home/aferrerferrer/Desktop/prueba"  #[PasspathVideoCreation]
    set pathVideo $path_Video 
    set postFileList $arg2
    set name_file $arg3
#[PassListFiles]
    GiD_Process Mescape Files ReadMultiple -rebuildIndex:0 $postFileList escape escape escape escape escape
    GiD_Process 'Zoom Frame Mescape escape escape escape escape
    GiD_Process Results ContourFill GAUSS_YOUNG escape escape escape escape
    GiD_Process Utilities Variables PostUpdateWindows no escape escape escape escape Mescape
    GiD_Process Results ContOptions SetMinOptions MinColor White Mescape 
    GiD_Process Results ContOptions SetMaxOptions MaxColor Black Mescape 
    set i_time_step 1
    foreach iStep [GiD_Info post get all_steps time] {
        GiD_Process Results AnalysisSel time $iStep Mescape Mescape
#	GiD_Process Results ContOptions SetMinOptions ResetValue Mescape 
#        GiD_Process Results ContOptions SetMaxOptions ResetValue Mescape 
	GiD_Process DisplayStyle bodybound Mescape
        GiD_Process Results ContourFill GAUSS_YOUNG GAUSS_YOUNG Mescape
 	GiD_Process Results ContOptions SetMinOptions SetValue 0.00999999977648258209 Mescape
	GiD_Process Results ContOptions SetMaxOptions SetValue 1 Mescape
#	Gid_Process Utilities Variables PostUpdateWindows yes escape escape escape escape
	GiD_Process MEscape 'Hardcopy Options ShowLegends No MEscape
	GiD_Process MEscape 'Hardcopy Options ShowAxes No MEscape
	GiD_Process MEscape 'Hardcopy Options PrintLogo No MEscape
	GiD_Process 'Hardcopy PNG "${pathVideo}-${i_time_step}_${name_file}.png"
	incr i_time_step
    }
    GiD_Process Utilities Variables PostUpdateWindows yes escape escape escape escape 
    GiD_Process Results ContOptions SetMinOptions MinColor Standard Mescape 
    GiD_Process Results ContOptions SetMaxOptions MaxColor Standard Mescape 

}
