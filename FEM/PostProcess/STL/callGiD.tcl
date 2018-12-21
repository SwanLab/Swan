set path "/home/alex/Desktop/tclFiles/"
set tclFile "CreateSurfaceSTL.tcl"
source $path$tclFile 
set output "$path/oe" 
set inputFile "/home/alex/git-repos/FEM-MAT-OO/Output/GrippingTriangleFine_Case_1_1_1/GrippingTriangleFine_Case_1_1_1_12.flavia.res"
set meshFile "$path/oe" 
set gidProjectName "$path/oe" 
set gidBasPath "/opt/GiDx64/13.0.2/templates/DXF.bas" 
CreateSurfaceSTL $inputFile $output $meshFile $gidProjectName $gidBasPath 
