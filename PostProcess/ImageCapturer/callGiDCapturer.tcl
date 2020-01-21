set path "/home/alex/git-repos/Swan/PostProcess/ImageCapturer/"
set tclFile "CaptureImage3.tcl"
source $path$tclFile 
set output /home/alex/git-repos/Swan/Output/RectangleMacroTriangleFineFine/TotalSmoothedRectangleMesh3Epsilon6 
set inputFile Output/RectangleMacroTriangleFineFine/RectangleMacroTriangleFineFine6.flavia.res
CaptureImage $inputFile $output 
