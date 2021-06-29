proc CaptureImage {pathFileVar outputFileVar} {
    set pathFile $pathFileVar
    set outputFile $outputFileVar
    GiD_Process PostProcess
    GiD_Process Mescape Files ReadMultiple -rebuildIndex:0 $pathFile escape escape escape escape escape
    GiD_Process 'Zoom Frame Mescape escape escape escape escape
    GiD_Process Results ContourFill LevelSet escape escape escape escape
    GiD_Process Results ContOptions ColorRamp DefineMap 3 #000000 #000000 #000000 Mescape
    GiD_Process Results ContOptions SetMaxOptions SetValue 0 Mescape Mescape
    GiD_Process Results ContOptions SetMaxOptions OutMaxColor White Mescape
    GiD_Process Results ContourFill LevelSet escape escape escape escape
	GiD_Process MEscape 'Hardcopy Options ShowLegends No MEscape
	GiD_Process MEscape 'Hardcopy Options ShowAxes No MEscape
	GiD_Process MEscape 'Hardcopy Options PrintLogo No MEscape
	GiD_Process 'Hardcopy PNG "$outputFile.png"
    GiD_Process Quit
}