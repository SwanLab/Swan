set path "/home/alex/git-repos/Swan/PostProcess/ImageCapturer/"
set tclFile "CaptureImage.tcl"
source $path$tclFile 
set output /home/alex/git-repos/Swan/Output/MySmoothCorner/MySmoothCorner1 
set inputFile Output/MySmoothCorner/MySmoothCorner1.flavia.res
CaptureImage $inputFile $output 
