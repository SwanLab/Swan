set path "/home/alex/git-repos/Swan/PostProcess/STL/"
set tclFile "GenerateMesh.tcl"
source $path$tclFile 
set gidProjectName "$path/sampleMesh" 
GenerateMesh $gidProjectName 
