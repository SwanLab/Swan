set path "/home/alex/git-repos/Swan/PostProcess/ImageCapturer/"
set tclFile "CaptureImage.tcl"
source $path$tclFile 
set output /home/alex/git-repos/Swan/Output/SmallCircleQ2/SmallCircleQ21 
set inputFile Output/SmallCircleQ2/SmallCircleQ21.flavia.res
CaptureImage $inputFile $output 
