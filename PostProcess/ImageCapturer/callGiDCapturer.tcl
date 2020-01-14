set path "/home/alex/git-repos/Swan/PostProcess/ImageCapturer/"
set tclFile "CaptureImage3.tcl"
source $path$tclFile 
set output /home/alex/git-repos/Swan/Output/SquareMacroTriangleFine/SmoothedCircleMesh2Epsilon6 
set inputFile Output/SquareMacroTriangleFine/SquareMacroTriangleFine6.flavia.res
CaptureImage $inputFile $output 
