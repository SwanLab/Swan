proc CaptureImage {pathFileVar outputFileVar} {
    set pathFile $pathFileVar
    set outputFile $outputFileVar
    GiD_Process PostProcess
    GiD_Process Mescape Files ReadMultiple -rebuildIndex:0 $pathFile escape escape escape escape escape
    GiD_Process 'Zoom Frame Mescape escape escape escape escape
    GiD_Process Mescape Results ContOptions NumberOfColor 50  Mescape     
    GiD_Process Results ContOptions ColorRamp Tangent Mescape
    GiD_Process escape escape escape escape escape Results ContourFill Density Density Mescape
    GiD_Process Results contoptions setmaxoptions setvalue 1 Mescape
    GiD_Process results contoptions setminoptions setvalue 0 Mescape
    GiD_Process Mescape DisplayStyle Body_Bound Mescape
    GiD_Process Results ContourFill Density Density Mescape
	GiD_Process MEscape 'Hardcopy Options ShowLegends No MEscape
	GiD_Process MEscape 'Hardcopy Options ShowAxes No MEscape
	GiD_Process MEscape 'Hardcopy Options PrintLogo No MEscape
	GiD_Process 'Hardcopy PNG "$outputFile.png"
    GiD_Process Quit
}