set path "/home/alex/git-repos/Swan/PostProcess/ImageCapturer/"
set tclFile "CaptureImage.tcl"
source $path$tclFile 
set output /home/alex/git-repos/Swan/Output//OptimalSuperEllipse0 
set inputFile Output/0.flavia.res
CaptureImage $inputFile $output 
