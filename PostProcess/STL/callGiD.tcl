set path "/home/joseantonio/Documentos/GitHub/Swan/PostProcess/STL/"
set tclFile "CreateSurfaceSTL.tcl"
source $path$tclFile 
set output "$path/sampleMesh" 
set inputFile "/home/joseantonio/Documentos/GitHub/Swan/samplemesh.flavia.res"
set meshFile "$path/sampleMesh" 
set gidProjectName "$path/sampleMesh" 
set gidBasPath "/home/joseantonio/GiDx64/gid-15.0.4/templates/DXF.bas" 
CreateSurfaceSTL $inputFile $output $meshFile $gidProjectName $gidBasPath 
