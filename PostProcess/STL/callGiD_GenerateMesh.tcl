set path "/home/gerard/Documents/GitHub/Swan/PostProcess/STL/"
set tclFile "GenerateMesh.tcl"
source $path$tclFile 
set gidProjectName "$path/sampleMesh" 
GenerateMesh $gidProjectName 
