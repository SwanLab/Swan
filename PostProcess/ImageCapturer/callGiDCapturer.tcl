set path "/home/alex/git-repos/Swan/PostProcess/ImageCapturer/"
set tclFile "CaptureImage3.tcl"
source $path$tclFile 
set output /home/alex/git-repos/MicroStructurePaper/OptimalMicroForStressCase3Phi7Print0 
set inputFile Output/OptimalMicroForStressCase3Phi7Print/OptimalMicroForStressCase3Phi7Print0.flavia.res
CaptureImage $inputFile $output 
