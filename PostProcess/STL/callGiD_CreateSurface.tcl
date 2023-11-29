set path "/home/alex/git-repos/Swan/PostProcess/STL/"
set tclFile "CreateSurfaceNew.tcl"
source $path$tclFile 
set output "$path/sampleMesh" 
set inputFile "/home/alex/git-repos/Swan/InnerMeshCreator_File.flavia.res"
set meshFile "$path/sampleMesh" 
set gidProjectName "$path/sampleMesh" 
set gidBasPath "/home/alex/bin/GiDx64/15.1.3d/templates/DXF.bas" 
CreateSurfaceNew $inputFile $output $meshFile $gidProjectName $gidBasPath 
