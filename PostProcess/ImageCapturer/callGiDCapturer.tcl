set path "/home/alex/git-repos/Swan/PostProcess/ImageCapturer/"
set tclFile "CaptureImageColor.tcl"
source $path$tclFile 
set output SmoothedCircleMesh1Epsilon1 
set inputFile Output/SquareMacroTriangle/SquareMacroTriangle1.flavia.res
CaptureImage $inputFile $output 
