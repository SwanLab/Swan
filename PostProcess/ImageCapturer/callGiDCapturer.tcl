set path "/home/alex/git-repos/Swan/PostProcess/ImageCapturer/"
set tclFile "CaptureImage3.tcl"
source $path$tclFile 
set output SmoothedCircleMesh2Epsilon1 
set inputFile Output/SquareMacroTriangleFine/SquareMacroTriangleFine1.flavia.res
CaptureImage $inputFile $output 
