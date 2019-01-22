set path "/home/alex/git-repos/FEM-MAT-OO/FEM/PostProcess/ImageCapturer/"
set tclFile "CaptureImage.tcl"
source $path$tclFile 
set output /home/alex/Dropbox/Amplificators/Images/RectangularInclusion0 
set inputFile AmplificatorsPdependency.mat
CaptureImage $inputFile $output 
