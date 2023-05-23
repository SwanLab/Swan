set path "/home/joseantonio/Documentos/GitHub/Swan/PostProcess/STL/"
set tclFile "ExportSTL.tcl"
source $path$tclFile 
set input "$path/sampleMesh.gid" 
set output "$path/sampleMeshFile.stl" 
set gidBasPath "/home/joseantonio/GiDx64/gid-15.0.4/templates/STL.bas" 
ExportSTL $input $output $gidBasPath 
