set path "/home/alex/Desktop/tclFiles/"
set tclFile "ExportSTL.tcl"
source $path$tclFile 
set input "$path/oe.gid" 
set output "$path/oeFile.stl" 
set gidBasPath "/opt/GiDx64/13.0.2/templates/STL.bas" 
ExportSTL $input $output $gidBasPath 
