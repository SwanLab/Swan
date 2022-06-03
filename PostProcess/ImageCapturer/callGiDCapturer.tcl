set path "/home/joseantonio/Documentos/GitHub/Swan/PostProcess/ImageCapturer"
set tclFile "CaptureImage3.tcl"
source $path$tclFile 
set output "testgidpic" 
set inputFile "test_anisotropy_cantilever"
CaptureImage $inputFile $output 
