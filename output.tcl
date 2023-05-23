set path "/home/joseantonio/Documentos/GitHub/Swan/PostProcess/STL/"
set tclFile "ExtrudeSTL_template.tcl"
source $path$tclFile 
set input "/home/joseantonio/Documentos/sampleMesh.gid" 
set output "$path/outputSTL.stl" 
set gidBasPath "/home/joseantonio/GiDx64/gid-15.0.4/templates/STL.bas" 
ExportSTL $input $output $gidBasPath 
