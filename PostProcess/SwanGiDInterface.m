classdef SwanGiDInterface < handle
    
    properties (Access = public)
        swanPath
        gidPath
        tclPath
    end
    
    properties (Access = private)
        
    end
    
    methods (Access = public)

        function obj = SwanGiDInterface(cParams)
            obj.swanPath = '/home/ton/Github/Swan/';
            obj.gidPath  = '/home/ton/GiDx64/gid-16.1.2d/';
            obj.tclPath  = [obj.swanPath, 'PostProcess/STL/'];
        end

        function extrudeAndExport(obj)
            resultsFile = '/home/ton/Github/Swan/hellouNou.flavia.res';
            obj.writeTclFile(resultsFile);
            obj.writeExportTclFile();
            command = [obj.gidPath,'gid_offscreen -offscreen -t "source ',obj.tclPath,'callGiD.tcl"'];
            system(command);
            command = [obj.gidPath,'gid_offscreen -offscreen -t "source ',obj.tclPath,'callGiD2.tcl"'];
            system(command);
        end

    end

    methods (Access = private)

        function writeTclFile(obj, resultsFile)
            tclFile = [obj.tclPath,'callGiD.tcl'];
            stlFileTocall = 'CreateSurfaceNew.tcl';
            gidBasPath = [obj.gidPath,'templates/DXF.bas'];
            fid = fopen(tclFile,'w+');
            fprintf(fid,['set path "',obj.tclPath,'"\n']);
            fprintf(fid,['set tclFile "',stlFileTocall,'"\n']);
            fprintf(fid,['source $path$tclFile \n']);
            fprintf(fid,['set output "$path/sampleMesh" \n']);
            fprintf(fid,['set inputFile "',resultsFile,'"\n']);
            fprintf(fid,['set meshFile "$path/sampleMesh" \n']);
            fprintf(fid,['set gidProjectName "$path/sampleMesh" \n']);
            fprintf(fid,['set gidBasPath "',gidBasPath,'" \n']);
            fprintf(fid,['CreateSurfaceNew $inputFile $output $meshFile $gidProjectName $gidBasPath \n']);
            fclose(fid);
        end

        function writeExportTclFile(obj)
            tclFile = [obj.tclPath,'.tcl'];
            stlFileTocall = 'ExportSTL.tcl';
            fid = fopen(tclFile,'w+');
            gidBasPath = [obj.gidPath,'templates/STL.bas'];
            fprintf(fid,['set path "',obj.tclPath,'"\n']);
            fprintf(fid,['set tclFile "',stlFileTocall,'"\n']);
            fprintf(fid,['source $path$tclFile \n']);
            fprintf(fid,['set input "$path/sampleMesh.gid" \n']);
            fprintf(fid,['set output "$path/sampleMeshFile.stl" \n']);
            fprintf(fid,['set gidBasPath "',gidBasPath,'" \n']);
            fprintf(fid,['ExportSTL $input $output $gidBasPath \n']);
            fclose(fid);
        end

    end

end

