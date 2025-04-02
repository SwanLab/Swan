proc CaptureImage {pathFileVar outputFileVar} {
    set pathFile $pathFileVar
    set outputFile $outputFileVar
    GiD_Process PostProcess
    GiD_Process Mescape Files ReadMultiple -rebuildIndex:0 $pathFile escape escape escape escape escape
    GiD_Process 'Zoom Frame Mescape escape escape escape escape
    GiD_Process Mescape Results ContOptions NumberOfColor 50  Mescape 
    GiD_Process Mescape DisplayStyle Body_Bound Mescape
    GiD_Process Results ContourFill DesignVar1 DesignVar1 Mescape
    GiD_Process Results contoptions setmaxoptions setvalue 1 Mescape
    GiD_Process Results contoptions setminoptions setvalue 0 Mescape
    GiD_Process Results ContOptions setminOptions MinColor White Mescape
    GiD_Process Results ContOptions setmaxOptions MaxColor Black Mescape
    GiD_Process Results ContOptions ColorRamp Tangent Mescape
    GiD_Process Utilities Redisplay escape escape    
	GiD_Process MEscape 'Hardcopy Options ShowLegends No MEscape
	GiD_Process MEscape 'Hardcopy Options ShowAxes No MEscape
	GiD_Process MEscape 'Hardcopy Options PrintLogo No MEscape
	GiD_Process 'Hardcopy PNG "$outputFile.png"
    GiD_Process Results ContOptions SetMinOptions MinColor Standard Mescape
    GiD_Process Results ContOptions SetMaxOptions MaxColor Standard Mescape
    GiD_Process Results ContOptions ColorRamp Tangent MEscape
#    GiD_Process Quit
}
