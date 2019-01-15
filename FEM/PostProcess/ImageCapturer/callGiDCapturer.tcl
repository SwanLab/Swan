set path "/home/alex/git-repos/FEM-MAT-OO/FEM/PostProcess/ImageCapturer/"
set tclFile "CaptureImage.tcl"
source $path$tclFile 
set output /home/alex/Dropbox/Amplificators/Images/RectangularInclusion0 
set inputFile /home/alex/git-repos/FEM-MAT-OO/Output/AmplificatorsPdependency.mat/AmplificatorsPdependency.mat0.flavia.res
CaptureImage $inputFile $output 
