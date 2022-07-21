set path "/home/alex/git-repos/Swan/PostProcess/ImageCapturer/"
set tclFile "CaptureImage3.tcl"
source $path$tclFile 
set output /home/alex/Dropbox/MaterialDesign/CC/Poisson/PoissonIter162 
set inputFile /media/alex/MyPassport/MaterialDesign/CStar/NegPoissonNoPerimeter5x6/HorizontalMaterialDesign162.flavia.res
CaptureImage $inputFile $output 
