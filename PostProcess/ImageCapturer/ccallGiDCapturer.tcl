set path "/home/alex/git-repos/Swan/PostProcess/ImageCapturer/c"
set tclFile "CaptureImage2.tcl"
source $path$tclFile 
set output /home/alex/git-repos/Swan/Output/VigergauzMicroStructure/VigergauzMicroStructure1 
set inputFile Output/VigergauzMicroStructure/VigergauzMicroStructure1.flavia.res
CaptureImage $inputFile $output 
