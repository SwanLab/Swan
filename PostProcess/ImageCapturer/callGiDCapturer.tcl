set path "/home/alex/git-repos/Swan/PostProcess/ImageCapturer/"
set tclFile "CaptureImage3.tcl"
source $path$tclFile 
set output /home/alex/git-repos/Swan/Output/VigdergauzPrinting/VigdergauzPrinting10 
set inputFile Output/VigdergauzPrinting/VigdergauzPrinting10.flavia.res
CaptureImage $inputFile $output 
