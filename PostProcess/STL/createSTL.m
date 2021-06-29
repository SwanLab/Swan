function createSTL
pathTcl = '/home/alex/Desktop/tclFiles/';
gidPath = '/home/alex/GiDx64/15.0.1/';
resultsFile = '/home/alex/git-repos/FEM-MAT-OO/Output/GrippingTriangleFine_Case_1_1_1/GrippingTriangleFine_Case_1_1_1_12.flavia.res';
writeTclFile(pathTcl,gidPath,resultsFile)
writeExportTclFile(pathTcl,gidPath)
command = [gidPath,'gid -t "source ',pathTcl,'callGiD.tcl"'];
unix(command);

delete('oe.msh')
delete('oe.png')
delete('oe.res')
delete('oe.vv')
delete('oe')

command = [gidPath,'gid -t "source callGiD2.tcl"'];
unix(command);
end

function writeExportTclFile(pathTcl,gidpath)
tclFile = 'callGiD2.tcl';
stlFileTocall = 'ExportSTL.tcl';
fid = fopen(tclFile,'w+');
gidBasPath = [gidpath,'templates/STL.bas'];
fprintf(fid,['set path "',pathTcl,'"\n']);
fprintf(fid,['set tclFile "',stlFileTocall,'"\n']);
fprintf(fid,['source $path$tclFile \n']);
fprintf(fid,['set input "$path/oe.gid" \n']);
fprintf(fid,['set output "$path/oeFile.stl" \n']);
fprintf(fid,['set gidBasPath "',gidBasPath,'" \n']);
fprintf(fid,['ExportSTL $input $output $gidBasPath \n']);
fclose(fid);
end

function writeTclFile(pathTcl,gidpath,resultsFile)
tclFile = 'callGiD.tcl';
stlFileTocall = 'CreateSurfaceSTL.tcl';
gidBasPath = [gidpath,'templates/DXF.bas'];
fid = fopen(tclFile,'w+');
fprintf(fid,['set path "',pathTcl,'"\n']);
fprintf(fid,['set tclFile "',stlFileTocall,'"\n']);
fprintf(fid,['source $path$tclFile \n']);
fprintf(fid,['set output "$path/oe" \n']);
fprintf(fid,['set inputFile "',resultsFile,'"\n']);
fprintf(fid,['set meshFile "$path/oe" \n']);
fprintf(fid,['set gidProjectName "$path/oe" \n']);
fprintf(fid,['set gidBasPath "',gidBasPath,'" \n']);
fprintf(fid,['CreateSurfaceSTL $inputFile $output $meshFile $gidProjectName $gidBasPath \n']);
fclose(fid);
end