set path "/home/gerard/Documents/GitHub/Swan/PostProcess/STL/"
set tclFile "CreateSurfaceNew.tcl"
source $path$tclFile 
set output "$path/sampleMesh" 
set inputFile "/home/gerard/Documents/GitHub/Swan/InnerMeshCreator_File.flavia.res"
set meshFile "$path/sampleMesh" 
set gidProjectName "$path/sampleMesh" 
set gidBasPath "/home/GiDx64/gid-16.0.6 templates/DXF.bas" 
CreateSurfaceNew $inputFile $output $meshFile $gidProjectName $gidBasPath 
