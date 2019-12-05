proc CaptureImage {pathFileVar outputFileVar} {
    set pathFile $pathFileVar
    set outputFile $outputFileVar
    GiD_Process PostProcess
    GiD_Process Mescape Files ReadMultiple -rebuildIndex:0 $pathFile escape escape escape escape escape
    GiD_Process 'Zoom Frame Mescape escape escape escape escape
    GiD_Process Results ContourFill LevelSet escape escape escape escape
    GiD_Process Utilities Variables PostUpdateWindows no escape escape escape escape Mescape
    GiD_Process Results ContOptions SetMinOptions MinColor White Mescape 
    GiD_Process Results ContOptions SetMaxOptions MaxColor Black Mescape 
	GiD_Process DisplayStyle bodybound Mescape
    GiD_Process Results ContourFill LevelSet LevelSet Mescape  
    GiD_Process Results ContOptions NumberOfColor 2  Mescape
    GiD_Process Results ContOptions SetMinOptions MinColor White Mescape
    GiD_Process Results ContOptions SetMaxOptions MaxColor White Mescape
    GiD_Process Results ContOptions SetMaxOptions SetValue 0 Mescape
    GiD_Process Results ContOptions ColorRamp Tangent Mescape
    GiD_Process Results ContOptions SetMinOptions OutMinColor Black Mescape
	GiD_Process MEscape 'Hardcopy Options ShowLegends No MEscape
	GiD_Process MEscape 'Hardcopy Options ShowAxes No MEscape
	GiD_Process MEscape 'Hardcopy Options PrintLogo No MEscape
	GiD_Process 'Hardcopy PNG "$outputFile.png"
    GiD_Process Quit
}
      
