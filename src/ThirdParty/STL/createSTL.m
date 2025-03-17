function createSTL
pathTcl = '/home/ton/Github/Swan/PostProcess/STL/';
gidPath = '/home/ton/GiDx64/gid-16.1.2d/';
resultsFile = '/home/ton/Github/Swan/Output/hellothere/hellothere1.flavia.res';

writeTclFile(pathTcl,gidPath,resultsFile)
writeExportTclFile(pathTcl,gidPath)
command = [gidPath,'gid_offscreen -offscreen -t "source ',pathTcl,'callGiD.tcl"'];
system(command);

command = [gidPath,'gid_offscreen -offscreen -t "source ',pathTcl,'callGiD2.tcl"'];
system(command);
end

function writeTclFile(pathTcl,gidpath,resultsFile)
tclFile = [pathTcl,'callGiD.tcl'];
stlFileTocall = 'CreateSurfaceSTL.tcl';
gidBasPath = [gidpath,'templates/DXF.bas'];
fid = fopen(tclFile,'w+');
fprintf(fid,['set path "',pathTcl,'"\n']);
fprintf(fid,['set tclFile "',stlFileTocall,'"\n']);
fprintf(fid,['source $path$tclFile \n']);
fprintf(fid,['set output "$path/sampleMesh" \n']);
fprintf(fid,['set inputFile "',resultsFile,'"\n']);
fprintf(fid,['set meshFile "$path/sampleMesh" \n']);
fprintf(fid,['set gidProjectName "$path/sampleMesh" \n']);
fprintf(fid,['set gidBasPath "',gidBasPath,'" \n']);
fprintf(fid,['CreateSurfaceSTL $inputFile $output $meshFile $gidProjectName $gidBasPath \n']);
fclose(fid);
end

function writeExportTclFile(pathTcl,gidpath)
tclFile = [pathTcl,'.tcl'];
stlFileTocall = 'ExportSTL.tcl';
fid = fopen(tclFile,'w+');
gidBasPath = [gidpath,'templates/STL.bas'];
fprintf(fid,['set path "',pathTcl,'"\n']);
fprintf(fid,['set tclFile "',stlFileTocall,'"\n']);
fprintf(fid,['source $path$tclFile \n']);
fprintf(fid,['set input "$path/sampleMesh.gid" \n']);
fprintf(fid,['set output "$path/sampleMeshFile.stl" \n']);
fprintf(fid,['set gidBasPath "',gidBasPath,'" \n']);
fprintf(fid,['ExportSTL $input $output $gidBasPath \n']);
fclose(fid);
end