set path "/home/alex/git-repos/Swan/PostProcess/ImageCapturer/"
set tclFile "CaptureImage2.tcl"
source $path$tclFile 
set output /home/alex/git-repos/Swan/Output/OptimalSuperEllipseAveragingSuperEllipseRho0_5Txi0_472Txi1Phi20/OptimalSuperEllipseAveragingSuperEllipseRho0_5Txi0_472Txi1Phi200 
set inputFile Output/OptimalSuperEllipseAveragingSuperEllipseRho0_5Txi0_472Txi1Phi20/OptimalSuperEllipseAveragingSuperEllipseRho0_5Txi0_472Txi1Phi200.flavia.res
CaptureImage $inputFile $output 
