set path "/home/ton/Github/Swan/PostProcess/STL/"
set tclFile "CreateSurfaceSTL.tcl"
source $path$tclFile 
set output "$path/sampleMesh" 
set inputFile "/home/ton/Github/Swan/Output/hellothere/hellothere1.flavia.res"
set meshFile "$path/sampleMesh" 
set gidProjectName "$path/sampleMesh" 
set gidBasPath "/home/ton/GiDx64/gid-16.1.2d/templates/DXF.bas" 
CreateSurfaceSTL $inputFile $output $meshFile $gidProjectName $gidBasPath 
