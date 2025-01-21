set path "/home/ton/Github/Swan/PostProcess/STL/"
set tclFile "ExportSTL.tcl"
source $path$tclFile 
set input "$path/sampleMesh.gid" 
set output "$path/sampleMeshFile.stl" 
set gidBasPath "/home/ton/GiDx64/gid-16.1.2d/templates/STL.bas" 
ExportSTL $input $output $gidBasPath 
