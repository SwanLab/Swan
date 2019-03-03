set path "/home/alex/git-repos/FEM-MAT-OO/PostProcess/ImageCapturer/"
set tclFile "CaptureImage.tcl"
source $path$tclFile 
set output /home/alex/Dropbox/Amplificators/Images/FourthOrderAmplificator0 
set inputFile Output/FourthOrderAmplificator/FourthOrderAmplificator0.flavia.res
CaptureImage $inputFile $output 
