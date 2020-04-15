set path "/home/alex/git-repos/Swan/PostProcess/ImageCapturer/"
set tclFile "CaptureImage2.tcl"
source $path$tclFile 
set output /home/alex/git-repos/Swan/Output/OptimalSuperEllipseExamplePapermax/OptimalSuperEllipseExamplePapermax0 
set inputFile Output/OptimalSuperEllipseExamplePapermax/OptimalSuperEllipseExamplePapermax0.flavia.res
CaptureImage $inputFile $output 
