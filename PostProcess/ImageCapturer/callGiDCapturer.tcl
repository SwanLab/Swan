set path "/home/alex/git-repos/Swan/PostProcess/ImageCapturer/"
set tclFile "CaptureImage2.tcl"
source $path$tclFile 
set output /home/alex/git-repos/Swan/Output/AveragingSuperEllipseRho0_5Txi0_472Txi1Phi2Q18_68/AveragingSuperEllipseRho0_5Txi0_472Txi1Phi2Q18_680 
set inputFile Output/AveragingSuperEllipseRho0_5Txi0_472Txi1Phi2Q18_68/AveragingSuperEllipseRho0_5Txi0_472Txi1Phi2Q18_680.flavia.res
CaptureImage $inputFile $output 
