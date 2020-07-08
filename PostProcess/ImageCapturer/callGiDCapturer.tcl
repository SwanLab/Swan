set path "/home/alex/git-repos/Swan/PostProcess/ImageCapturer/"
set tclFile "CaptureImage3.tcl"
source $path$tclFile 
set output /home/alex/git-repos/MicroStructurePaper/ExamplePaper0_1Print0 
set inputFile Output/ExamplePaper0_1Print/ExamplePaper0_1Print0.flavia.res
CaptureImage $inputFile $output 
