set path "/home/alex/git-repos/Swan/PostProcess/ImageCapturer/"
set tclFile "CaptureImage3.tcl"
source $path$tclFile 
set output /home/alex/git-repos/Swan/Output/SquareMacroTriangle/SmoothedCircleMesh1Epsilon5 
set inputFile Output/SquareMacroTriangle/SquareMacroTriangle5.flavia.res
CaptureImage $inputFile $output 
