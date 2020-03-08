set path "/home/alex/git-repos/Swan/PostProcess/ImageCapturer/"
set tclFile "CaptureImage2.tcl"
source $path$tclFile 
set output /home/alex/git-repos/Swan/Output/OptimalSuperEllipse5/OptimalSuperEllipse55 
set inputFile Output/OptimalSuperEllipse5/OptimalSuperEllipse55.flavia.res
CaptureImage $inputFile $output 
