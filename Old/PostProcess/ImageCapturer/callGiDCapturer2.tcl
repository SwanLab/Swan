set path "/home/alex/git-repos/FEM-MAT-OO/FEM/PostProcess/ImageCapturer/"
set tclFile "CaptureImage.tcl"
source $path$tclFile 
set output "/home/alex/Desktop/FotoMicro" 
set inputFile "/home/alex/git-repos/FEM-MAT-OO/Output/AmplificatorTensorForInclusionDensity/AmplificatorTensorForInclusionDensity_2.flavia.res"
CaptureImage $inputFile $output
