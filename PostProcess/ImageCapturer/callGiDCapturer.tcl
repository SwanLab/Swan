set path "/home/alex/git-repos/Swan/PostProcess/ImageCapturer/"
set tclFile "CaptureImage3.tcl"
source $path$tclFile 
set output /home/alex/git-repos/Swan/Output/CircleMacro/SmoothedCircleMesh1Epsilon1 
set inputFile Output/CircleMacro/CircleMacro1.flavia.res
CaptureImage $inputFile $output 
