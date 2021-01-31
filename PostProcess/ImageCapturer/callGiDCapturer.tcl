set path "/home/alex/git-repos/Swan/PostProcess/ImageCapturer/"
set tclFile "CaptureImageColor.tcl"
source $path$tclFile 
set output RelativeSmoothedRectangleMeshColor1Epsilon6 
set inputFile Output/RectangleMacroTriangle/RectangleMacroTriangle6.flavia.res
CaptureImage $inputFile $output 
